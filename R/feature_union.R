#' Union of Selected Features from Multiple Feature Selection Results
#'
#' Extracts and combines all \code{feature_selection} vectors found within one or more
#' nested list objects (such as outputs from \code{FIAs}, \code{FIAs_W},
#' \code{Cluster_priority_lasso}, \code{boruta_regression},
#' \code{auto_abess}, \code{glmnet_regression}, or \code{mmpc_regression}).
#' Returns the unique union of all selected feature names.
#'
#' Additionally, evaluates the predictive performance of the union feature set using
#' k-fold cross-validated Random Forest regression (\code{ranger}) on the provided dataset.
#'
#' @param ... One or more list objects. Each list may contain nested elements with
#'   \code{feature_selection} entries, which are expected to be character vectors of feature names.
#' @param data A data frame with the first column as the target variable and the remaining columns as features.
#' @param metrics Character vector of mlr3 regression performance measures to evaluate (default: "regr.rsq").
#' @param seed Integer random seed for reproducibility of the model training.
#' @param n_folds Number of folds for cross-validation (default is 10).
#'
#' @return A named list containing:
#' \describe{
#'   \item{\code{union_scores}}{Data frame with iteration-wise cross-validation scores for the union feature set.}
#'   \item{\code{union_aggregates}}{Aggregated performance metrics across all folds for the union feature set.}
#'   \item{\code{union_features}}{Returns a vector consisting of the union of features constructed from all inputs.}
#'   \item{\code{feature_ratio}}{Ratio of selected features in the union relative to all features in \code{data}.}
#'   \item{\code{Prior_informed_matrix}}{Binary matrix indicating feature selection status by each input list (rows = methods, columns = features).}
#' }
#'
#' @examples
#' \dontrun{
#' #You can select any combination of results to construct a feature union.
#' #Uses some results to construct a union (excluding FIAs).
#' union1 <- feature_union(FIAs_W_result$results$FIAs_W,abess_result,Boruta_result,PL_result,EN_result,mmpc_result,data = WW1768)
#' #Uses all results to construct a union.
#' union2 <- feature_union(FIAs_W_result,abess_result,Boruta_result,PL_result,EN_result,mmpc_result,data = WW1768)
#' }
#'
#' @export
feature_union <- function(..., data, metrics = "regr.rsq", seed = 1, n_folds = 10) {
  `%>%` <- magrittr::`%>%`

  find_feature_selection <- function(lst) {
    result <- list()
    if (is.list(lst)) {
      if ("feature_selection" %in% names(lst)) {
        fr <- lst[["feature_selection"]]
        if (is.vector(fr)) {
          result <- c(result, list(fr))
        }
      }
      for (element in lst) {
        if (is.list(element)) {
          result <- c(result, find_feature_selection(element))
        }
      }
    }
    return(result)
  }

  all_lists <- list(...)
  all_features <- unlist(lapply(all_lists, function(lst) {
    unlist(find_feature_selection(lst))
  }))
  all_features <- unique(all_features)

  target_var <- names(data)[1]
  feature_ratio <- length(all_features) / (length(data)-1)

  union_task <- data %>%
    dplyr::select(all_of(c(target_var, all_features)))
  task <- mlr3::as_task_regr(union_task, target = target_var)

  base::set.seed(seed)
  learner <- mlr3::lrn("regr.ranger")
  resampling <- mlr3::rsmp("cv", folds = n_folds)
  rr <- mlr3::resample(task, learner, resampling, store_models = FALSE)

  union_scores <- rr$score(mlr3::msrs(metrics)) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(algorithm = "union") %>%
    dplyr::select(c("algorithm", "iteration", "learner_id"), dplyr::any_of(c("rsq", metrics)))

  union_aggregates <- rr$aggregate(mlr3::msrs(metrics)) %>%
    as.list() %>%
    tibble::as_tibble() %>%
    dplyr::mutate(algorithm = "union") %>%
    dplyr::select(c("algorithm", dplyr::any_of(c("rsq", metrics))))

  init_matrix <- do.call(rbind, lapply(all_lists, function(lst) {
    selected <- unlist(find_feature_selection(lst))
    as.integer(all_features %in% selected)
  }))
  rownames(init_matrix) <- paste0("method_", seq_len(length(all_lists)))
  colnames(init_matrix) <- all_features

  return(list(
    union_scores = union_scores,
    union_aggregates = union_aggregates,
    union_features = all_features,
    feature_ratio = feature_ratio,
    Prior_Guided_matrix = init_matrix
  ))
}
