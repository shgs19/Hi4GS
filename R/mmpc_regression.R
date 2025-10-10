#' Feature Selection Using MMPC regression
#'
#' Applies the MMPC (Max-Min Parents and Children) algorithm for feature selection in a regression context.
#' For each specified value of \code{max_k}, it selects features, builds a regression model using \code{ranger},
#' and evaluates performance through cross-validation
#'
#' @param data Data frame format, the first column is the target variable (trait value), and the remaining columns are SNP features. The encoding format of these features can be 0, 1, 2 or -1, 0, 1.
#' @param seed Numeric format, random seed for random forest model training.
#' @param n_folds Numeric format, default is 10, specifies the cross validation fold of the random forest model.
#' @param max_k_range Integer vector. A sequence of values for the \code{max_k} parameter in the MMPC algorithm. Each value defines the maximum size of the conditioning set used in conditional independence tests. Increasing \code{max_k} allows more complex dependency structures to be detected but increases computation time.
#' @param threshold Numeric. The significance level (usually between 0 and 1, e.g., 0.05) used in the MMPC algorithm to determine whether a variable should be included. A lower value makes the selection more stringent.
#' @param ncores Integer. Number of processor cores to use for parallel computation during feature selection. If \code{NULL}, it defaults to one less than the number of available physical cores.
#' @param metrics Character vector specifying the performance metrics used to evaluate the random forest model.
#' The final variable selection is based on "regr.rsq".
#' Defaults to 21 common regression metrics, including "regr.rsq", "regr.mse", "regr.rmse", among others.
#' @param output_file Character string, default "MMPC_result.xlsx". Name of the Excel file to save the MMPC feature statistics and model evaluation results.
#'
#' @details
#' The MMPC (Max-Min Parents and Children) algorithm is a constraint-based feature selection method
#' that identifies relevant variables by iteratively testing conditional independencies,
#' aiming to find the direct causes and effects (parents and children) of a target variable.
#' The parameter \code{max_k} controls the maximum size of the conditioning set in these tests,
#' allowing the detection of increasingly complex dependency structures as \code{max_k} increases,
#' at the cost of greater computational time.
#'
#' This function applies MMPC for feature selection across a range of \code{max_k} values.
#' For each \code{max_k}, it selects relevant features from the input dataset,
#' then fits a random forest regression model using \code{ranger} on the selected features.
#' The model's predictive performance is evaluated by k-fold cross-validation according to multiple metrics.
#'
#' @return A list containing:
#' \describe{
#'   \item{scores}{A data frame of cross-validation scores for each support size and fold.}
#'   \item{aggregates}{A data frame of aggregated scores for each support size.}
#'   \item{feature_selection}{A character vector of confirmed important features selected by MMPC.}
#'   \item{best_max_k}{The \code{max_k} value yielding the best R-squared.}
#'   \item{rsq}{The best cross-validated R-squared value.}
#'   \item{time_mins}{A numeric value indicating the total runtime of the function in minutes.}
#'   \item{memory_bytes}{A numeric value representing the approximate memory usage (in bytes) during execution, computed via \code{pryr::mem_used()}.}
#' }
#'
#' @examples
#' \dontrun{
#' mmpc_result <- mmpc_regression(data = WW1768, n_folds = 10,max_k_range = 1:10,threshold = 0.05,seed = 1)
#' }
#'
#' @references
#' Tsamardinos I, Brown LE, Aliferis CF (2006) The max-min hill-climbing Bayesian network structure learning algorithm. Machine learning 65(1):31-78.
#' \url{https://doi.org/10.1007/s10994-006-6889-7}
#'
#' @seealso \code{\link{boruta_regression},\link{Cluster_priority_lasso},\link{auto_abess},\link{auto_abess}}}
#'
#' @export
mmpc_regression <- function(
    data,
    seed = 1,
    n_folds,
    max_k_range = 1:10,
    threshold = 0.05,
    ncores = NULL,
    metrics = c("regr.rsq", "regr.mse", "time_both", "regr.rmse", "regr.bias",
                "regr.ktau", "regr.mae", "regr.mape", "regr.maxae", "regr.medae",
                "regr.medse", "regr.msle", "regr.pbias", "regr.rae", "regr.rmsle",
                "regr.rrse", "regr.rse", "regr.sae", "regr.smape", "regr.srho", "regr.sse"),
    output_file = "MMPC_result.xlsx"
) {
  `%>%` <- magrittr::`%>%`
  start_time <- Sys.time()
  gc()
  start_mem <- pryr::mem_used()

  if (is.null(ncores)) {
    ncores <- parallel::detectCores(logical = FALSE)
    if (is.na(ncores) || ncores < 1) {
      ncores <- 1
    } else {
      ncores <- max(1, ncores - 1)
    }
  }

  target = names(data)[1]
  y <- data[[target]]
  x <- data %>% dplyr::select(-all_of(target))

  mmpc_list <- purrr::map(max_k_range, function(k) {
    message("\nRunning max_k = ", k)
    mmpc_model <- MXM::MMPC(
      target = y,
      dataset = x,
      ncores = ncores,
      max_k = k,
      threshold = threshold
    )

    selected_vars <- names(x)[mmpc_model@selectedVars]
    selected_data <- data %>%
      dplyr::select(all_of(c(target, selected_vars)))

    mmpc_tk <- mlr3::as_task_regr(selected_data, target = target)

    set.seed(seed)
    learner <- mlr3::lrn("regr.ranger")
    resampling <- mlr3::rsmp("cv", folds = n_folds)
    rr <- resample(mmpc_tk, learner, resampling, store_models = TRUE)

    mmpc_scores <- rr$score(mlr3::msrs(metrics)) %>%
      tibble::as_tibble() %>%
      dplyr::mutate(algorithm = "mmpc", max_k = k) %>%
      dplyr::select(c("algorithm", "max_k", "iteration", "learner_id"), dplyr::any_of(c("rsq", metrics)))

    mmpc_aggregates <- rr$aggregate(mlr3::msrs(metrics)) %>%
      as.list() %>%
      tibble::as_tibble() %>%
      dplyr::mutate(algorithm = "mmpc", max_k = k) %>%
      dplyr::select(c("algorithm", "max_k"), dplyr::any_of(c("rsq", metrics)))

    list(scores = mmpc_scores, aggregates = mmpc_aggregates, selected_vars = selected_vars)
  }) %>% purrr::compact()

  end_time <- Sys.time()
  gc()
  end_mem <- pryr::mem_used()
  total_time <- round(as.numeric(difftime(end_time, start_time, units = "mins")), 3)
  total_bytes <- end_mem - start_mem

  final_scores <- purrr::map_dfr(mmpc_list, "scores")
  final_aggregates <- purrr::map_dfr(mmpc_list, "aggregates")

  best_idx <- which.max(final_aggregates$rsq)
  rsq <- final_aggregates$rsq[best_idx]
  best_max_k <- final_aggregates$max_k[best_idx]
  feature_selection <- mmpc_list[[best_idx]]$selected_vars

  writexl::write_xlsx(
    list(
      "mmpc Scores" = as.data.frame(final_scores),
      "mmpc Aggregates" = as.data.frame(final_aggregates),
      "time_memory" = tibble::tibble(
        rsq = rsq,time_mins = total_time,memory_bytes = total_bytes)),
       path = output_file)

  list(
    scores = final_scores,
    aggregates = final_aggregates,
    feature_selection = feature_selection,
    best_max_k = best_max_k,
    rsq = rsq,
    time_mins = total_time,
    memory_bytes = total_bytes
  )
}
