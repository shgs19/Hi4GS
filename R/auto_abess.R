#' Automated Feature Selection and Model Evaluation Using abess
#'
#' This function performs feature selection with the abess algorithm over a range of support sizes,
#' evaluates models via cross-validation, and identifies the best feature subset based on performance metrics.
#'
#' @param data Data frame format, the first column is the target variable (trait value), and the remaining columns are SNP features. The encoding format of these features can be 0, 1, 2 or -1, 0, 1.
#' @param start_k Integer. The initial number of features (support size) to try in the abess algorithm.
#' @param step Integer. The increment of support size for each iteration.
#' @param max_tries Integer. Maximum number of different support sizes (iterations) to attempt before stopping.
#' @param seed Numeric format, random seed for random forest model training.
#' @param n_folds Numeric format, default is 10, specifies the cross validation fold of the random forest model.
#' @param num.threads Integer. Number of threads to use for abess computation. If NULL, the function will automatically detect the number of available cores on the machine (using \code{parallel::detectCores(logical = TRUE)}) and use one less than the total number of cores.
#' @param max.splicing.iter Integer. Maximum number of splicing iterations in abess.
#' @param verbose Logical. Whether to print progress messages. Default is TRUE.
#' @param metrics Character vector specifying the performance metrics used to evaluate the random forest model.
#' The final variable selection is based on "regr.rsq".
#' Defaults to 21 common regression metrics, including "regr.rsq", "regr.mse", "regr.rmse", among others.
#' @param output_file Character string, default "abess_result.xlsx". Name of the Excel file to save the abess feature statistics and model evaluation results.
#'
#' @details
#' ABESS (Adaptive BEst Subset Selection) is a polynomial-time algorithm designed to efficiently
#' solve the best subset selection problem in high-dimensional regression. Unlike traditional
#' stepwise or Lasso-based methods, ABESS aims to find the optimal subset of features with strong
#' theoretical guarantees on selection accuracy and computational efficiency.
#'
#' This function automates the feature selection process by iteratively applying ABESS with
#' increasing support sizes (the number of selected features). At each iteration, it fits an ABESS
#' model with the specified \code{support.size} (equal to \code{current_k}), and evaluates the selected
#' feature subset by training a random forest regression model with cross-validation using
#' \code{n_folds} folds.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{parameters}{A list of parameters used in the function call.}
#'   \item{scores}{A data frame of cross-validation scores for each support size and fold.}
#'   \item{aggregates}{A data frame of aggregated scores for each support size.}
#'   \item{model_best}{A list containing the best model's scores, aggregate metrics, and selected features.}
#'   \item{rsq}{The best cross-validated R-squared value.}
#'   \item{time_mins}{Total execution time in minutes.}
#'   \item{memory_bytes}{Memory used during the execution.}
#'   \item{termination_reason}{Text explanation of why the loop stopped (limit reached, error, etc.).}
#' }
#'
#' @references
#' Zhu J, Wen C, Zhu J, et al (2020) A polynomial algorithm for best-subset selection problem. Proceedings of the National Academy of Sciences 117(52):33117-33123.
#' \url{https://doi.org/10.1073/pnas.2014241117}
#'
#' @examples
#' \dontrun{
#' # starting with 1 feature, increasing by 1 each iteration, up to 5 attempts, with verbose output enabled.
#' abess_result <- auto_abess(data = WW1768,start_k = 1,
#' step = 1,max_tries = 5,verbose = TRUE)
#' }
#'
#' @seealso \code{\link{boruta_regression},\link{Cluster_priority_lasso},\link{glmnet_regression},\link{mmpc_regression}}}
#'
#' @export
auto_abess <- function(
    data,
    start_k = 1,
    step = 1,
    max_tries = 1000,
    seed = 1,
    n_folds = 10,
    num.threads = NULL,
    max.splicing.iter = 50,
    verbose = TRUE,
    metrics = c("regr.rsq", "regr.mse", "time_both", "regr.rmse", "regr.bias",
                "regr.ktau", "regr.mae", "regr.mape", "regr.maxae", "regr.medae",
                "regr.medse", "regr.msle", "regr.pbias", "regr.rae", "regr.rmsle",
                "regr.rrse", "regr.rse", "regr.sae", "regr.smape", "regr.srho", "regr.sse"),
    output_file = "abess_result.xlsx"
) {
  `%>%` <- magrittr::`%>%`
  start_time <- Sys.time(); gc(); start_mem <- pryr::mem_used()

  target = names(data)[1]
  if (ncol(data) < 2) stop("Data must contain at least 1 feature column")
  if (!target %in% names(data)) stop("Target variable does not exist in the data")
  if (is.null(num.threads)) {
    num.threads <- parallel::detectCores(logical = TRUE)-1
  }

  n_features <- ncol(data) - 1
  results <- list()
  current_k <- start_k
  continue <- TRUE

  for (attempt in 1:max_tries) {
    if (!continue || current_k > n_features) break
    if (verbose) message("\nTesting support.size = ", current_k)
    result <- tryCatch({
      emb_filter <- mlr3filters::flt("selected_features",
                        learner = mlr3::lrn("regr.abess",
                                      support.size = current_k,
                                      c.max = current_k,
                                      max.splicing.iter = max.splicing.iter,
                                      nfolds = n_folds,
                                      warm.start = TRUE,
                                      num.threads = num.threads))

      task <- mlr3::as_task_regr(data, target = target)
      emb_filter$calculate(task)
      features <- data.table::as.data.table(emb_filter)$feature[1:current_k]

      if (length(features) < current_k) stop("Actual selected features count (", length(features), ") < requested (", current_k, ")")
      if (any(is.na(features))) stop("NA features detected")

      abess_tk <- data %>% dplyr::select(all_of(c(target, features)))

      base::set.seed(seed)
      learner <- mlr3::lrn("regr.ranger")
      resampling <-  mlr3::rsmp("cv", folds = n_folds)
      rr <- mlr3::resample(task = mlr3::as_task_regr(abess_tk, target = target),
                     learner = learner,
                     resampling = resampling,
                     store_models = FALSE)

      list(
        score = rr$score(msrs(metrics)) %>%
          tibble::as_tibble() %>%
          dplyr::mutate(support_size = current_k) %>%
          dplyr::relocate(support_size, .before = iteration),
        aggregate = rr$aggregate(msrs(metrics)),
        support_size = current_k,
        feature_selection = features
      )
    }, error = function(e) {
      if (verbose) message("! Testing interrupted: ", e$message)
      continue <<- FALSE
      return(NULL)
    })

    if (!is.null(result)) {
      results[[as.character(current_k)]] <- result
      current_k <- current_k + step
    }
  }

  end_time <- Sys.time(); gc(); end_mem <- pryr::mem_used()
  total_time <- round(as.numeric(difftime(end_time, start_time, units = "mins")), 3)
  total_bytes <- end_mem - start_mem

  abess_aggregate <- dplyr::bind_rows(lapply(results, function(x) {
      tibble::as_tibble_row(x$aggregate) %>%
        dplyr::mutate(support_size = x$support_size)
    }))

  model_best <- if (nrow(abess_aggregate) > 0) {
    if ("rsq" %in% names(abess_aggregate)) {
      idx <- which.max(abess_aggregate[["rsq"]])
      results[[idx]]
    } else {
      NULL
    }
  }

  scores_df = if (length(results) > 0) dplyr::bind_rows(lapply(results, `[[`, "score")) else tibble()
  abess_score = scores_df %>% dplyr::select(c(task_id, support_size, iteration, learner_id),any_of(c("rsq",metrics)))
  rsq = model_best[["aggregate"]][['rsq']]

  writexl::write_xlsx(
    list(
      "abess Scores" = as.data.frame(abess_score),
      "abess Aggregates" = as.data.frame(abess_aggregate),
      "best Scores" = as.data.frame(model_best[["score"]]),
      "best Aggregates" = as.data.frame(as.list(model_best[["aggregate"]])),
      "time_memory" = tibble::tibble(
        time_mins = total_time,
        memory_bytes = total_bytes)),path = output_file)

  list(
    parameters = list(
      target = target,
      start_k = start_k,
      step = step,
      seed = seed,
      n_folds = n_folds,
      max_tries = max_tries,
      num.threads = num.threads,
      max.splicing.iter = max.splicing.iter,
      verbose = verbose,
      metrics = metrics
    ),
    scores = abess_score,
    aggregates = abess_aggregate,
    model_best = model_best,
    rsq = rsq,
    time_mins = total_time,
    memory_bytes = total_bytes,
    termination_reason = ifelse(current_k > n_features,
     paste("Reached feature limit:", n_features),
     ifelse(continue, "Completed all attempts", "Terminated due to error"))
  )
}
