#' Feature Selection Using Elastic Net
#'
#' Performs feature selection using cross-validated elastic net regression (`cv_glmnet`) and evaluates the predictive performance
#' of the selected features using a random forest model (`ranger`) with k-fold cross-validation.
#' The function returns both per-fold performance scores and aggregated performance metrics.
#'
#' @param data Data frame format, the first column is the target variable (trait value), and the remaining columns are SNP features. The encoding format of these features can be 0, 1, 2 or -1, 0, 1.
#' @param seed Numeric format, random seed for random forest model training.
#' @param n_folds Numeric format, default is 10, specifies the cross validation fold of the random forest model.
#' @param metrics Character vector specifying the performance metrics used to evaluate the random forest model.
#' The final variable selection is based on "regr.rsq".
#' Defaults to 21 common regression metrics, including "regr.rsq", "regr.mse", "regr.rmse", among others.
#' @param output_file Character string, default "glmnet_result.xlsx". Name of the Excel file to save the elastic net feature statistics and model evaluation results.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{scores}{A data frame of cross-validation scores for each support size and fold.}
#'   \item{aggregates}{A data frame of aggregated scores for each support size.}
#'   \item{feature_selection}{A character vector of confirmed important features selected by elastic net regression.}
#'   \item{rsq}{The best cross-validated R-squared value.}
#'   \item{time_mins}{A numeric value indicating the total runtime of the function in minutes.}
#'   \item{memory_bytes}{A numeric value representing the approximate memory usage (in bytes) during execution, computed via \code{pryr::mem_used()}.}
#' }
#'
#' @details
#' The elastic net algorithm is a regularized regression technique that integrates L1 (Lasso) and L2 (Ridge) penalties to simultaneously perform variable selection and coefficient shrinkage. By tuning the mixing parameter alpha, elastic net effectively handles correlated predictors and selects groups of related features, addressing limitations of using Lasso or Ridge alone. The optimal regularization parameter lambda is chosen through cross-validation to minimize prediction error.
#' This function first applies cross-validated elastic net regression (using the `cv_glmnet` filter from `mlr3`) to score and identify informative SNP features predictive of the target trait. Selected features are then used to build a random forest regression model (`ranger`), which undergoes k-fold cross-validation to evaluate predictive performance with multiple metrics such as R-squared, MSE, RMSE, and bias.
#' The process also tracks execution time and memory consumption. Results—including per-fold scores, aggregated metrics, and the final set of selected features—are saved to an Excel file to facilitate result interpretation and reproducibility. By combining elastic net’s sparse, robust feature selection with random forest’s nonlinear modeling capability, this workflow provides an efficient and reliable framework for SNP selection and evaluation in genomic prediction.
#'
#' @examples
#' \dontrun{
#' EN_result <- glmnet_regression(data = WW1768,seed = 1,n_folds = 10)
#' }
#'
#' @references
#' Friedman JH, Hastie T, Tibshirani R (2010) Regularization paths for generalized linear models via coordinate descent. Journal of statistical software 33(1):1-22.
#' \url{https://doi.org/10.18637/jss.v033.i01}
#'
#' Zou H, Hastie T (2005) Regularization and variable selection via the elastic net. Journal of the Royal Statistical Society Series B: Statistical Methodology 67(2):301-320.
#' \url{https://doi.org/10.1111/j.1467-9868.2005.00527.x}
#'
#' @seealso \code{\link{boruta_regression},\link{Cluster_priority_lasso},\link{auto_abess},\link{mmpc_regression}}}
#'
#' @export
glmnet_regression <- function(
    data,
    seed = 1,
    n_folds = 10,
    metrics = c("regr.rsq", "regr.mse", "time_both", "regr.rmse", "regr.bias",
                "regr.ktau", "regr.mae", "regr.mape", "regr.maxae", "regr.medae",
                "regr.medse", "regr.msle", "regr.pbias", "regr.rae", "regr.rmsle",
                "regr.rrse", "regr.rse", "regr.sae", "regr.smape", "regr.srho", "regr.sse"),
    output_file = "glmnet_result.xlsx"
) {
  `%>%` <- magrittr::`%>%`
  start_time <- Sys.time(); gc(); start_mem <- pryr::mem_used()
  target = names(data)[1]
  score_threshold = 1

  if (!target %in% names(data)) stop("Target variable not found in data")
  if (ncol(data) < 2) stop("Data must contain at least 1 feature column")

  task <- as_task_regr(data, target = target)

  flt <- mlr3filters::flt("selected_features",
             learner = mlr3::lrn("regr.cv_glmnet",
                           nfolds = n_folds,
                           s = "lambda.min"))
  flt$calculate(task)
  score_vec <- flt$scores
  selected_features <- tibble::tibble(
    feature = names(score_vec),
    score = as.numeric(score_vec)
  ) %>%
    dplyr::filter(score >= score_threshold) %>%
    dplyr::arrange(desc(score))

  if (nrow(selected_features) == 0) {
    warning("No features selected with score >= ", score_threshold)
    return(NULL)
  }

  selected_data <- data %>%
    dplyr::select(all_of(c(target, selected_features$feature)))

  glmnet_tk <- mlr3::as_task_regr(selected_data, target = target)

  set.seed(seed)
  learner <- mlr3::lrn("regr.ranger")
  resampling <- mlr3::rsmp("cv", folds = n_folds)
  rr <- resample(glmnet_tk, learner, resampling, store_models = FALSE)

  end_time <- Sys.time(); gc(); end_mem <- pryr::mem_used()
  total_time <- round(as.numeric(difftime(end_time, start_time, units = "mins")), 3)
  total_bytes <- end_mem - start_mem

  glmnet_scores <- rr$score(mlr3::msrs(metrics)) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(algorithm = "glmnet") %>%
    dplyr::select(c(algorithm, iteration, learner_id), any_of(c("rsq",metrics)))

  glmnet_aggregates <- rr$aggregate(mlr3::msrs(metrics)) %>%
    as.list() %>%
    tibble::as_tibble() %>%
    dplyr::mutate(algorithm = "glmnet") %>%
    dplyr::select(algorithm, any_of(c("rsq",metrics)))

  rsq = glmnet_aggregates$rsq

  writexl::write_xlsx(
    list("glmnet Scores" = as.data.frame(glmnet_scores),
    "glmnet Aggregates" = as.data.frame(glmnet_aggregates),
    "time_memory" = tibble::tibble(
    rsq = rsq,time_mins = total_time,memory_bytes = total_bytes)),
    path = output_file)


  list(
    scores = glmnet_scores,
    aggregates = glmnet_aggregates,
    feature_selection = selected_features$feature,
    rsq = rsq,
    time_mins = total_time,
    memory_bytes = total_bytes
  )
}
