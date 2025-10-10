#' Genomic Selection with BWGS and mlr3 Frameworks
#'
#' This function performs genomic selection (GS) by integrating classical statistical methods available in the \code{BWGS} package
#' (e.g., GBLUP, Bayesian methods) with flexible machine-learning regression learners from the \code{mlr3} ecosystem.
#' It runs k-fold cross-validation, computes prediction accuracy via Pearson’s correlation coefficient (PCC), and saves results.
#'
#' @param data A data frame where the first column is the target trait (phenotypic response) and the rest are SNP marker values.
#' Row names correspond to genotype identifiers. If absent or purely sequential,row names will be auto-generated as “variety1”, “variety2”, etc.
#' @param methods Character vector of prediction methods to evaluate. Can include BWGS-supported methods (e.g. "GBLUP","EGBLUP","BA","BB","BC", etc.)
#' and any regression learner key from \code{mlr3extralearners} (e.g., "regr.ranger").
#' Defaults to \code{c("GBLUP","EGBLUP","BA","BB","BC")}.
#' @param seed Integer. Random seed for reproducible fold splits and model training.
#' @param n_folds Integer. Number of folds for cross-validation (default: 10).
#' @param output_file File path to save results in Excel format (.xlsx).Default: "gs_results.xlsx".
#'
#' @details
#' The function first separates the genotype matrix (SNPs) and phenotype vector (trait values).
#' Genotypes are optionally imputed (via \code{MNI} within BWGS) for BWGS-supported methods.
#' The dataset is split into \code{n_folds} folds for cross-validation.
#'
#' For each method and fold:
#' \enumerate{
#'   \item The model is trained on the training set (genotype + phenotype).
#'   \item Predictions are generated for the test set.
#'   \item Pearson correlation (PCC) between observed and predicted trait values is computed as the performance metric.
#' }
#'
#' The function supports both BWGS native prediction methods and machine learning models
#' accessible through the \code{mlr3} framework (including learners from \code{mlr3extralearners}).
#' For BWGS methods, genotype imputation and prediction are handled internally,
#' while \code{mlr3} learners are trained and predicted using the \code{TaskRegr} interface.
#'
#' Results are aggregated across folds for each method, including:
#' \itemize{
#'   \item Fold-level PCC values,
#'   \item Per-genotype predictions (observed vs. predicted),
#'   \item Average PCC across folds per method.
#' }
#'
#' @return an R data frame (nested table) with the following structure:
#' \describe{
#' \item{\code{methods}}{Prediction methods used during evaluation}
#' \item{\code{avg_pcc}}{Average PCC value for each method using K-fold cross validation}
#' \item{\code{cv_folds}}{A nested data frame of fold-level PCC results. }
#' \item{\code{raw_data}}{A nested data frame of observed and predicted values for each method in each fold. }
#' }
#'
#' @references
#' Lang M, Binder M, Richter J, et al (2019) mlr3: A modern object-oriented machine learning framework in R. Journal of Open Source Software 4(44):1903.
#' \url{https://doi.org/10.21105/joss.01903}
#'
#' Bischl, Bernd, et al (2024) Applied machine learning using mlr3 in R. CRC Press.
#' \href{https://books.google.co.uk/books?hl=zh-CN&lr=&id=5wrsEAAAQBAJ&oi=fnd&pg=PP1&dq=Bischl,+B.,+et+al.+(2024).+Applied+Machine+Learning+Using+mlr3+in+R.&ots=8cVZQFDnfp&sig=lxMJLXgOP4vt8Q97n23OxKPbC04&redir_esc=y#v=onepage&q=Bischl%2C%20B.%2C%20et%20al.%20(2024).%20Applied%20Machine%20Learning%20Using%20mlr3%20in%20R.&f=false}{books.google}
#'
#' Charmet G, Tran LG, Auzanneau J, et al (2020) BWGS: AR package for genomic selection and its application to a wheat breeding programme. PLoS One 15(4):e0222733.
#' \url{https://doi.org/10.1371/journal.pone.0232422}
#'
#' @examples
#' \dontrun{
#' # Run Genomic Selection with BWGS and mlr3
#' Hi4GS.data <- WW1768 %>% select(1,PGWO_result$feature_selection)
#' gs_results <- gs_bwgs_mlr3(
#'   data = Hi4GS.data,
#'   methods = c("GBLUP","regr.ranger"),
#'   seed = 1,
#'   n_folds = 10,
#'   output_file = "gs_results.xlsx"
#' )
#'
#' # Access average PCC per method
#' gs_results$avg_pcc
#'
#' # Inspect predictions for a specific method
#' gs_results$raw_data[[1]]
#' }
#'
#' @export
gs_bwgs_mlr3 <- function(data, methods = c("GBLUP","EGBLUP","BA", "BB", "BC"), seed = 1, n_folds = 10, output_file = "gs_results.xlsx") {

data <- as.data.frame(data)
Trait <- colnames(data)[1]

if (is.null(rownames(data)) || !all(rownames(data) == as.character(seq_len(nrow(data))))) {
  geno <- data[, -1] %>% as.matrix()
  pheno <- data[[1]]
  names(pheno) <- rownames(data)
} else {
  message("Row names are missing. Automatically generating: variety1, variety2, ...")
  rownames(data) <- paste0('variety', 1:nrow(data))
  geno <- data[, -1] %>% as.matrix()
  pheno <- data[[1]]
  rownames(geno) <- rownames(data)
  names(pheno) <- rownames(data)
}

n <- length(pheno)
set.seed(seed)
folds <- sample(rep(1:n_folds, length.out = n))

bwgs_methods <- c("GBLUP", "EGBLUP", "RR", "LASSO", "EN","BRR", "BL", "BA", "BB", "BC","RKHS", "RF", "SVM", "BRNN")
mlr3_methods <- mlr3extralearners::lrns()$keys()[grepl("^regr\\.", mlr3extralearners::lrns()$keys())]
all_methods <- methods

bwgs_df <- purrr::map_dfr(all_methods, function(m) {
  fold_results <- purrr::map_dfr(1:n_folds, function(k) {
    train_idx <- which(folds != k)
    test_idx  <- which(folds == k)
    geno_train <- geno[train_idx, ]
    pheno_train <- pheno[train_idx]
    geno_test <- geno[test_idx, ]

    if (m %in% bwgs_methods) {
      pred <- BWGS::bwgs.predict(
        geno_train = geno_train,
        pheno_train = pheno_train,
        geno_target = geno_test,
        MAP = "NULL",
        predict.method = m
      ) %>% as.data.frame()
      y_pred <- pred[[1]]
    } else {
      task <- mlr3::TaskRegr$new(id = "gs", backend = data.frame(y = pheno_train, geno_train), target = "y")
      learner <- mlr3::lrn(m)
      learner$train(task)
      test_data <- data.frame(y = pheno[test_idx], geno_test)
      task_test <- mlr3::TaskRegr$new(id = "gs_test", backend = test_data, target = "y")
      y_pred <- learner$predict(task_test)$response
    }

    dplyr::tibble(
      variety = names(pheno)[test_idx],
      trait = Trait,
      observed = pheno[test_idx],
      predicted = y_pred,
      methods = m,
      folds = k
    )
  })
}) %>% dplyr::arrange(folds, methods)

cv_folds <- bwgs_df %>%
  dplyr::group_by(methods, trait, folds) %>%
  dplyr::summarise(pcc = stats::cor(observed, predicted, use = "pairwise.complete.obs", method = "pearson"),
                   .groups = "drop")

main_table <- cv_folds %>%
  dplyr::group_by(methods, trait) %>%
  dplyr::summarise(avg_pcc = mean(pcc, na.rm = TRUE),.groups = "drop")

gs_nest <- main_table %>%
  dplyr::left_join(
    bwgs_df %>%
      dplyr::mutate(algorithm = methods) %>%
      tidyr::nest(raw_data = c(variety, trait, observed, predicted, algorithm, folds)),
    by = c("methods")
  ) %>%
  dplyr::left_join(
    cv_folds %>%
      dplyr::mutate(algorithm = methods) %>%
      tidyr::nest(cv_folds = c(algorithm, trait, folds, pcc)),
    by = c("methods")
  )

writexl::write_xlsx(
  list(
    bwgs_df = bwgs_df,
    cv_folds = cv_folds,
    main_table = main_table
  ),
  output_file
)

message("\nExcel file saved to: ", normalizePath(output_file))

return(gs_nest)
}
