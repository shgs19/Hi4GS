#' Feature Importance Algorithms with Incremental Feature Selection
#'
#' A comprehensive SNP screening approach was implemented by integrating multiple feature importance measures, including mRMR, SNP's PVE, Pearson correlation between SNPs and traits, permutation-based random forest importance, absolute coefficients from ridge regression, and GWAS p-values, in combination with incremental feature selection.
#'
#' @param data Data frame format, the first column is the target variable (trait value), and the remaining columns are SNP features. The encoding format of these features can be 0, 1, 2 or -1, 0, 1.
#' @param algorithms Character vector format, default is c("mrmr", "PVE", "Cor","RF","ridge","gwas"), case insensitive, these algorithms perform feature selection and input the obtained features into RF to return the evaluation index. In particular, when gwas is selected, it is necessary to read the two data frame parameters gwas.myY and gwas.myGM related to the algorithm.
#' @param n_features Numeric format, with the default set to the total number of features. This parameter defines the maximum number of features to be selected, effectively determining the value of K in the Top K features. The features are ranked by importance and sequentially added to the model one by one until the number of selected features reaches the specified n_features (i.e., incremental feature selection). This process will use the maximum rsq of the model performance in 1:n_features as the criterion for feature selection.
#' @param n_folds Numeric format, default is 10, specifies the cross validation fold of the random forest model.
#' @param output_path Output result file, default is the path where getwd() is located.
#' @param seed Numeric format, random seed for random forest model training.
#' @param mrmr.threads Numeric format, the number of threads of the mRMR algorithm, the default is the maximum number of threads on the computer minus one.
#' @param ridge.lambda A continuous integer vector specifying the search range for lambda in the ridge algorithm. The default value is 1:300.
#' @param ridge.metrics Character format, performance evaluation index of ridge algorithm, default is "regr.rsq".
#' @param ridge.evals Numeric format, the maximum number of evaluations during hyperparameter tuning, that is, the maximum number of lambda values the tuner will test. If NULL, the number of evaluations defaults to the range (max - min) of \code{ridge.lambda}.
#' @param gwas.myY A data frame required when using the GWAS algorithm. The first column should contain genotype names, and the second column should contain the trait values.
#' The genotype names in the first column must match the genotype names in the feature column of \code{data}.
#' @param gwas.myGM A data frame required when using the GWAS algorithm. It must have three columns: the first for SNP identifiers, the second for chromosome numbers, and the third for the corresponding physical positions.
#' @param gwas.nPCA Numerical format, the number of principal components (PCA) included in the GWAS model, as used in the GAPIT package.
#' @param gwas.model Character format, one or more GWAS models supported by the GAPIT package, such as "FarmCPU", "MLM", "GLM", etc. Default is "FarmCPU".
#' @param metrics Character vector specifying the performance metrics used to evaluate the random forest model.
#' The final variable selection is based on "regr.rsq".
#' Defaults to 21 common regression metrics, including "regr.rsq", "regr.mse", "regr.rmse", among others.
#'
#' @details
#' This function performs feature selection by combining multiple feature importance algorithms and incremental feature selection (IFS).
#' First, features are ranked by importance using the specified algorithms.
#' Next, features are incrementally added into the predictive model (random forest by default) in order of their importance.
#' The model performance (evaluated via rsq) is tracked across 1:n_features, and the feature subset achieving the maximum rsq is selected as the final set of informative SNPs.
#'
#' @return A list with three components:
#' \describe{
#'   \item{data}{The input data frame.}
#'   \item{parameters}{A list of input parameters and settings.}
#'   \item{results}{A list of results for each specified feature importance algorithm, including the ranked features, model object, runtime, memory usage, and cross-validation outcomes.}
#' }
#'
#' @references
#' Bermingham ML, Pong-Wong R, Spiliopoulou A, et al (2015) Application of high-dimensional feature selection: evaluation for genomic prediction in man. Scientific reports 5(1):10312.
#' \url{https://doi.org/10.1038/srep10312}
#'
#' Wang J, Zhang Z (2021) GAPIT version 3: boosting power and accuracy for genomic association and prediction. Genomics, proteomics & bioinformatics 19(4):629-640.
#' \url{https://doi.org/10.1016/j.gpb.2021.08.005}
#'
#'  Lang M, Binder M, Richter J, et al (2019) mlr3: A modern object-oriented machine learning framework in R. Journal of Open Source Software 4(44):1903.
#' \url{https://doi.org/10.21105/joss.01903}
#'
#' Peng H, Long F, Ding C (2005) Feature selection based on mutual information criteria of max-dependency, max-relevance, and min-redundancy. IEEE Transactions on pattern analysis and machine intelligence 27(8):1226-1238.
#' \url{https://doi.org/10.1109/TPAMI.2005.159}
#'
#' Breiman L (2001) Random forests. Machine learning 45(1):5-32.
#' \url{https://doi.org/10.1023/A:1010933404324}
#'
#' Hoerl AE, Kennard RW (1970) Ridge regression: Biased estimation for nonorthogonal problems. Technometrics 12(1):55-67.
#' \url{https://doi.org/10.1080/00401706.1970.10488634}
#'
#' Friedman JH, Hastie T, Tibshirani R (2010) Regularization paths for generalized linear models via coordinate descent. Journal of statistical software, 33(1):1-22.
#' \url{https://doi.org/10.18637/jss.v033.i01}
#'
#' Wright MN, Wager S, Probst P (2020) Ranger: A fast implementation of random forests. R package version 0.12 1:730.
#' \url{https://doi.org/10.18637/jss.v077.i01}
#'
#' @examples
#' \dontrun{
#' # Running the Feature Importance Algorithms with Incremental Feature Selection
#' FIAs_result <- FIAs(
#'   data = WW1768,
#'   gwas.myY = gwas.myY,
#'   gwas.myGM = gwas.myGM,
#'   n_folds = 10,
#'   n_features = 500,
#'   ridge.lambda = 1:300,
#'   algorithms = c("mrmr", "PVE", "Cor", "RF", "ridge", "gwas")
#' )
#' rank_fias(data = FIAs_result)
#'}
#'
#' @seealso \code{\link{FIAs_W},\link{rank_fias}}
#'
#' @export
FIAs <- function(
    data,
    algorithms = c("mrmr", "PVE", "Cor","RF","ridge","gwas"),
    n_features = NULL,
    n_folds = 10,
    output_path = "./",
    seed = 1,
    mrmr.threads = NULL,
    ridge.lambda = 1:300,
    ridge.metrics = "regr.rsq",
    ridge.evals = NULL,
    gwas.myY,
    gwas.myGM,
    gwas.nPCA = 3,
    gwas.model = c("FarmCPU"),
    metrics = c("regr.rsq", "regr.mse", "time_both", "regr.rmse", "regr.bias",
                "regr.ktau", "regr.mae", "regr.mape", "regr.maxae", "regr.medae",
                "regr.medse", "regr.msle", "regr.pbias", "regr.rae", "regr.rmsle",
                "regr.rrse", "regr.rse", "regr.sae", "regr.smape", "regr.srho", "regr.sse")
) {
  `%>%` <- magrittr::`%>%`
  target_col <- names(data)[1]
  if (!is.data.frame(data)) stop("data must be a data.frame")
  if (ncol(data) < 2) stop("Data must contain at least 1 feature and 1 target column")
  if (!all(tolower(algorithms) %in% c("mrmr", "pve","cor","rf","ridge","gwas"))) stop("Supported algorithms: mrmr, PVE")
  if (!is.numeric(data[[target_col]])) stop("Target column must be numeric for regression")
  if (is.null(n_features)) {
    n_features <- ncol(data) - 1
    message("\nUsing all available features: n_features = ", n_features)
  }
  if (is.null(ridge.evals)) {
    ridge.evals <- diff(range(ridge.lambda)) + 1
  }
  if (is.null(mrmr.threads)) {
    mrmr.threads <- parallel::detectCores(logical = TRUE)-1
  }

  final_results <- list(
    data = data,
    parameters = list(
      n_features = n_features,
      n_folds = n_folds,
      algorithms = algorithms,
      target_variable = target_col,
      n_folds = n_folds,
      seed = seed,
      metrics = metrics,
      timestamp = Sys.time()
    ),
    results = list()
  )

  for (algo in algorithms) {
    message("\n[=== Processing ", toupper(algo), " Algorithm ===]")
    start_time <- Sys.time(); gc(); start_mem <- pryr::mem_used()

    algo_result <- switch(tolower(algo),
                          "mrmr" = {
                            message("\n[1/3] Calculating MRMR features...")
                            task <- mlr3::as_task_regr(data, target = target_col)
                            flt <- mlr3filters::flt("mrmr",threads = mrmr.threads)
                            flt$calculate(task)
                            feature_ranking <- data.table::as.data.table(flt)$feature

                            list(
                              model = flt$clone(),
                              mrmr_result = data.table::as.data.table(flt),
                              selected_features = feature_ranking
                            )
                          },

                          "cor" = {
                            message("\n[1/3] Calculating correlation features...")
                            task <- mlr3::as_task_regr(data, target = target_col)
                            flt <- mlr3filters::flt("correlation")
                            flt$calculate(task)
                            feature_ranking <- data.table::as.data.table(flt)$feature

                            list(
                              model = flt$clone(),
                              Cor_result = data.table::as.data.table(flt),
                              selected_features = feature_ranking
                            )
                          },

                          "rf" = {
                            message("\n[1/3] Calculating RF features...")
                            task <- mlr3::as_task_regr(data, target = target_col)
                            learner <- mlr3::lrn("regr.ranger", importance = "permutation")
                            learner$train(task)
                            importance <- learner$model$variable.importance

                            RF <- data.frame(
                              feature = names(importance),
                              importance = as.numeric(importance)
                            )
                            feature_ranking <- RF %>% dplyr::arrange(-importance) %>% dplyr::pull(feature)

                            list(
                              model = learner$clone(),
                              RF_result = as.data.table(RF),
                              selected_features = feature_ranking
                            )
                          },

                          "pve" = {
                            message("\n[1/3] Calculating PVE features...")
                            pca <- stats::prcomp(data[, -1, drop = FALSE])
                            pcs <- pca$x[, 1:3]
                            colnames(pcs) <- paste0("PC", 1:3)
                            SNP_genotypes <- cbind(as.data.frame(data[, -1]), pcs)
                            Trait <- data[[1]]

                            model2 <- stats::lm(Trait ~ PC1 + PC2 + PC3, data = SNP_genotypes)
                            r_squared2 <- summary(model2)$r.squared

                            PVE_result <- purrr::map_dfr(
                              seq_len(ncol(SNP_genotypes) - 3),
                              function(i) {
                                model1 <- lm(Trait ~ SNP_genotypes[[i]] + PC1 + PC2 + PC3, data = SNP_genotypes)
                                summ <- summary(model1)
                                effect <- summ$coefficients[2, "Estimate"]
                                p_value <- summ$coefficients[2, "Pr(>|t|)"]
                                r_squared1 <- summ$r.squared
                                PVE <- r_squared1 - r_squared2

                                tibble::tibble(
                                  SNP = colnames(SNP_genotypes)[i],
                                  Effect = effect,
                                  P_value = p_value,
                                  R2.snp = r_squared1,
                                  R2.pc = r_squared2,
                                  PVE = PVE
                                )
                              }
                            )

                            feature_ranking <- PVE_result %>%
                              dplyr::arrange(desc(PVE)) %>%
                              dplyr::pull(SNP)

                            list(
                              pca_model = pca,
                              pve_result = `rownames<-`(PVE_result, NULL),
                              selected_features = feature_ranking
                            )
                          },

                          "ridge" = {
                            message("\n[1/3] Calculating Ridge features...")

                            target_col <- names(data)[1]
                            task <- mlr3::as_task_regr(data, target = target_col)
                            learner <- mlr3::lrn("regr.glmnet",
                                                 alpha = 0,
                                                 standardize = TRUE)

                            search_space <- paradox::ps(
                              lambda = paradox::p_int(
                                lower = min(ridge.lambda),
                                upper = max(ridge.lambda)
                              )
                            )

                            tuner <- mlr3tuning::tnr("grid_search", resolution = ridge.evals)
                            instance <- mlr3tuning::TuningInstanceSingleCrit$new(
                              task = task,
                              learner = learner,
                              resampling = mlr3::rsmp("cv", folds = 5),
                              measure = mlr3::msr(ridge.metrics),
                              search_space = search_space,
                              terminator = mlr3tuning::trm("evals", n_evals = ridge.evals)
                            )

                            message("\n*Ridge* Tuning lambda hyperparameter...")
                            future::plan("multisession")
                            tuner$optimize(instance)

                            message("\n*Ridge* Training final model...")
                            learner$param_set$values <- instance$result_learner_param_vals
                            learner$train(task)

                            message("\n*Ridge* Calculating feature importance...")
                            coef_matrix <- glmnet::coef.glmnet(
                              learner$model,
                              s = instance$result$lambda
                            )

                            feature_importance <- data.frame(
                              feature = rownames(coef_matrix)[-1],
                              coefficient = as.vector(coef_matrix[-1, 1]),
                              importance = abs(as.vector(coef_matrix[-1, 1]))
                            ) %>% dplyr::arrange(desc(importance))
                            future::plan("sequential")

                            list(
                              model = learner,
                              Ridge_result = feature_importance,
                              selected_features = feature_importance$feature,
                              lambda = instance$result$lambda
                            )
                          },

                          "gwas" = {
                            message("[1/3] Calculating GWAS features...")
                            myY <- gwas.myY %>% data.frame()
                            myGD <- dplyr::bind_cols(myY[,1],data[,-1]) %>% data.frame()
                            names(myGD)[1] <- names(myY)[1]
                            myGM <- gwas.myGM %>% data.frame()
                            names(myGM) <- c("SNP","Chromosome","Position")

                            geno_matrix <- myGD[, -1, drop = FALSE]
                            unique_values <- unique(unlist(geno_matrix))
                            unique_values <- unique_values[!is.na(unique_values)]

                            if (all(unique_values %in% c(-1, 0, 1))) {
                              message("Genotype codes of -1/0/1 are detected and automatically converted to 0/1/2")
                              myGD <- myGD %>%
                                dplyr::mutate(across(-1, ~ .x + 1L))
                            } else if (all(unique_values %in% 0:2)) {
                              message("The genotype is in the 0/1/2 format")
                            } else {
                              invalid_values <- setdiff(unique_values, c(-1:2))
                              stop(
                                "Invalid genotype values detected: ",
                                paste(sort(invalid_values), collapse = ", "), "\n",
                                "Allowed genotype encoding formats: \n",
                                "1. Original encoding: -1(AA), 0(AC), 1(CC) \n",
                                "2. Standard encoding: 0(AA), 1(AC), 2(CC)\n",
                                "Input data must strictly follow one of these formats"
                              )
                            }

                            gwas_dir <- file.path(output_path, "GWAS_Results")
                            if (!dir.exists(gwas_dir)) {
                              message("Creating GWAS output directory: ", gwas_dir)
                              dir.create(gwas_dir, recursive = TRUE, showWarnings = FALSE)
                            }

                            myGAPIT <- withr::with_dir(gwas_dir, GAPIT::GAPIT(
                              Y = myY,
                              GD = myGD,
                              GM = myGM,
                              PCA.total = gwas.nPCA,
                              model = gwas.model,
                              Multiple_analysis = TRUE
                            ))

                            feature_ranking <- myGAPIT$GWAS %>%
                              dplyr::arrange(P.value) %>%
                              dplyr::pull(SNP)

                            list(
                              gwas_GAPIT = myGAPIT,
                              gwas_result = myGAPIT$GWAS,
                              selected_features = feature_ranking
                            )
                          }

    )

    actual_features <- length(algo_result$selected_features)
    algo_n_features <- min(n_features, actual_features)
    message("Adjusted n_features to ", algo_n_features,
            " (available features: ", actual_features, ")")

    message("[2/3] Running cross-validation...")
    pb <- txtProgressBar(min = 0, max = algo_n_features, style = 3)

    score_list <- vector("list", algo_n_features)
    aggregate_list <- vector("list", algo_n_features)

    for (count in seq_len(algo_n_features)) {
      setTxtProgressBar(pb, count)

      selected_cols <- c(target_col, algo_result$selected_features[1:count])
      selected_data <- data[, selected_cols]

      task_sub <- tryCatch({
        mlr3::as_task_regr(selected_data, target = target_col)
      }, error = function(e) {
        stop("Task creation failed at count ", count, ": ", e$message)
      })

      rr <- tryCatch({
        set.seed(seed)
        mlr3::resample(
          task = task_sub,
          learner = mlr3::lrn("regr.ranger"),
          resampling = mlr3::rsmp("cv", folds = n_folds),
          store_models = FALSE
        )
      }, error = function(e) {
        warning("Resampling failed at count ", count, ": ", e$message)
        NULL
      })

      if (!is.null(rr)) {

        scores <- rr$score(mlr3::msrs(metrics)) %>%
          tibble::as_tibble() %>%
          dplyr::mutate(gradient = count, algorithm = algo)

        aggregates <- rr$aggregate(mlr3::msrs(metrics)) %>%
          as.list() %>%
          tibble::as_tibble() %>%
          dplyr::mutate(gradient = count, algorithm = algo)

        score_list[[count]] <- scores
        aggregate_list[[count]] <- aggregates
      }
    }
    close(pb)

    message("[3/3] Finalizing ", algo, " results...")

    algo_scores <- data.table::rbindlist(score_list, fill = TRUE) %>%
      dplyr::select(c(algorithm, gradient, iteration, learner_id), any_of(c("rsq",metrics)))

    algo_aggregates <- data.table::rbindlist(aggregate_list, fill = TRUE) %>%
      dplyr::select(c(algorithm, gradient), any_of(c("rsq",metrics)))

    rsq_col <- grep("rsq", names(algo_aggregates), value = TRUE)
    optimal_k <- algo_aggregates[which.max(algo_aggregates[[rsq_col]]), ][["gradient"]]
    feature_selection <- head(algo_result$selected_features, optimal_k)

    end_time <- Sys.time(); gc(); end_mem <- pryr::mem_used()
    total_time <- round(as.numeric(difftime(end_time, start_time, units = "mins")), 3)
    total_bytes <- end_mem - start_mem

    final_results$results[[algo]] <- list(
      model = algo_result[!names(algo_result) %in% c("selected_features")],
      scores = algo_scores,
      aggregates = algo_aggregates,
      feature_ranking = algo_result$selected_features,
      feature_selection = feature_selection,
      time_mins = total_time,
      memory_bytes = total_bytes
    )

    output_file <- file.path(output_path, paste0(algo, "_results.xlsx"))
    writexl::write_xlsx(
      list(
        scores = as.data.frame(algo_scores),
        aggregates = as.data.frame(algo_aggregates),
        "time_memory" = tibble::tibble(
          time_mins = total_time,
          memory_bytes = total_bytes)), path = output_file)
    message("Results saved to:\n", normalizePath(output_file))
  }
  message("\n========== All Algorithms Completed ==========")
  invisible(final_results)
}
