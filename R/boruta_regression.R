#' Feature Selection Using Boruta Regression
#'
#' Perform feature selection based on Boruta and then test the performance in random forest.
#'
#' @param data Data frame format, the first column is the target variable (trait value), and the remaining columns are SNP features. The encoding format of these features can be 0, 1, 2 or -1, 0, 1.
#' @param maxRuns Maximum number of runs for the Boruta algorithm.
#' @param seed Numeric format, random seed for random forest model training.
#' @param n_folds Numeric format, default is 10, specifies the cross validation fold of the random forest model.
#' @param metrics Character vector specifying the performance metrics used to evaluate the random forest model.
#' The final variable selection is based on "regr.rsq".
#' Defaults to 21 common regression metrics, including "regr.rsq", "regr.mse", "regr.rmse", among others.
#' @param output_file Character string, default "Boruta_result.xlsx". Name of the Excel file to save the Boruta feature statistics and model evaluation results.
#'
#' @return A named list containing the following elements:
#' \describe{
#'   \item{boruta_model}{An object of class \code{"Boruta"}, the fitted Boruta model.}
#'   \item{boruta_stats}{A data frame of feature importance statistics returned by \code{attStats()}, including columns such as \code{meanImp}, \code{normHits}, and \code{decision}.}
#'   \item{scores}{A tibble of resampling performance scores for each fold, based on the selected metrics.}
#'   \item{aggregates}{A tibble of aggregated cross-validation metrics (e.g., mean R², RMSE) across all folds.}
#'   \item{feature_selection}{A character vector of confirmed important features selected by Boruta.}
#'   \item{time_mins}{A numeric value indicating the total runtime of the function in minutes.}
#'   \item{memory_bytes}{A numeric value representing the approximate memory usage (in bytes) during execution, computed via \code{pryr::mem_used()}.}
#' }
#'
#' @details
#' Boruta is a wrapper-based feature selection algorithm built around random forests.
#' It iteratively compares the importance of original features with importance achievable
#' at random (shadow features) to robustly identify truly relevant variables.
#'
#' This function performs feature selection using Boruta with a maximum number of iterations
#' specified by \code{maxRuns}. After Boruta confirms important features, a random forest
#' regression model is trained on the selected subset using cross-validation to evaluate
#' prediction performance.
#'
#' @examples
#' \dontrun{
#' # Run Boruta-based feature selection and evaluate with 10-fold CV
#' Boruta_result <- boruta_regression(data = WW1768, n_folds = 10, maxRuns = 200, seed = 1)
#' }
#'
#' @references
#' Kursa MB, Rudnicki WR (2010) Feature selection with the Boruta package. Journal of statistical software 36:1-13.
#' \url{https://doi.org/10.18637/jss.v036.i11}
#'
#' @seealso \code{\link{auto_abess},\link{Cluster_priority_lasso},\link{glmnet_regression},\link{mmpc_regression}}}
#'
#' @export
boruta_regression <- function(data, maxRuns, seed = 1, n_folds = 10,
                              metrics = c("regr.rsq", "regr.mse", "time_both", "regr.rmse", "regr.bias",
                                          "regr.ktau", "regr.mae", "regr.mape", "regr.maxae", "regr.medae",
                                          "regr.medse", "regr.msle", "regr.pbias", "regr.rae", "regr.rmsle",
                                          "regr.rrse", "regr.rse", "regr.sae", "regr.smape", "regr.srho", "regr.sse"),
                              output_file = "Boruta_result.xlsx") {
  `%>%` <- magrittr::`%>%`
  start_time <- Sys.time()
  gc()
  start_mem <- pryr::mem_used()

  target_var <- names(data)[1]
  boruta_model <- Boruta::Boruta(
    formula = stats::as.formula(paste(target_var, "~ .")),
    data = data,
    maxRuns = maxRuns
  )

  feature_stats <- Boruta::attStats(boruta_model)
  confirmed_features <- rownames(feature_stats[feature_stats$decision == "Confirmed", ])

  if (length(confirmed_features) == 0) {
    warning("No confirmed important features were found, Please adjust the maxRuns parameter")
    return(NULL)
  }

  Boruta_task <- data %>%
    dplyr::select(all_of(c(target_var, confirmed_features)))
  task <- mlr3::as_task_regr(Boruta_task, target = target_var)

  base::set.seed(seed)
  learner <- mlr3::lrn("regr.ranger")
  resampling <- mlr3::rsmp("cv", folds = n_folds)

  rr <- mlr3::resample(task, learner, resampling, store_models = FALSE)
  measures <- mlr3::msrs(metrics)

  end_time <- Sys.time()
  gc()
  end_mem <- pryr::mem_used()
  total_time <- round(as.numeric(difftime(end_time, start_time, units = "mins")), 3)
  total_bytes <- end_mem - start_mem

  Boruta_scores <- rr$score(mlr3::msrs(metrics)) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(algorithm = "Boruta") %>%
    dplyr::select(c(algorithm, iteration, learner_id), any_of(c("rsq",metrics)))

  Boruta_aggregates <- rr$aggregate(mlr3::msrs(metrics)) %>%
    as.list() %>%
    tibble::as_tibble() %>%
    dplyr::mutate(algorithm = "Boruta") %>%
    dplyr::select(algorithm, any_of(c("rsq",metrics)))

  feature_stats_rownames <- tibble::rownames_to_column(feature_stats, var = "Feature")
  writexl::write_xlsx(
    list(
      "Boruta stats" = feature_stats_rownames,
      "Boruta Scores" = as.data.frame(Boruta_scores),
      "Boruta Aggregates" = as.data.frame(Boruta_aggregates),
      "time_memory" = tibble::tibble(
        time_mins = total_time,memory_bytes = total_bytes)),path = output_file)

  list(
    boruta_model = boruta_model,
    boruta_stats = feature_stats,
    scores = Boruta_scores,
    aggregates = Boruta_aggregates,
    feature_selection = confirmed_features,
    time_mins = total_time,
    memory_bytes = total_bytes
  )
}
