#' Clustering-based Priority LASSO
#'
#' The block parameters of Priority LASSO are divided by three clustering methods: clara, kmeans, and gmm to improve model performance.
#'
#' @param data Data frame format, the first column is the target variable (trait value), and the remaining columns are SNP features. The encoding format of these features can be 0, 1, 2 or -1, 0, 1.
#' @param cluster_methods A character vector specifying the clustering methods to use. Options include `"clara"`, `"kmeans"`, and `"gmm"`. Default is all three.
#' @param k_values Optional numeric vector. Specifies the number(s) of clusters to test for each method. If `NULL`, defaults to `1:10` (or `2:10` for GMM) limited by the number of features.
#' @param clara.samples integer, say N, the number of samples to be drawn from the dataset. The default, N = 5, is rather small for historical (and now back compatibility) reasons and we recommend to set samples an order of magnitude larger.
#' @param seed Numeric format, random seed for random forest model training.
#' @param n_folds Numeric format, default is 10, specifies the cross validation fold of the random forest model.
#' @param output_file File name ending with .xlsx, output weight results, optimal weight, cross-validation performance, and feature ranking.
#'
#' @details
#' Priority LASSO is a feature selection method designed for situations
#' where features can be divided into meaningful blocks with a predefined priority order.
#' It performs block-wise LASSO penalization: features in higher-priority blocks
#' are entered into the model first, while lower-priority blocks are penalized more heavily.
#' This approach is especially useful in high-dimensional settings such as genomics,
#' where the number of features (\emph{p}) greatly exceeds the number of samples (\emph{n}).
#'
#' In its standard form, Priority LASSO requires the user to define the block structure of features.
#' However, when such a priori block information is unavailable,
#' clustering algorithms can be used to construct these blocks automatically.
#'
#' This function integrates clustering with Priority LASSO as follows:
#' \enumerate{
#'   \item Extract the SNP feature matrix from the input data (excluding the target variable).
#'   \item Apply one or more clustering algorithms (\code{"clara"}, \code{"kmeans"}, \code{"gmm"}) to group features into \emph{k} clusters (blocks).
#'   \item Use the resulting clusters as the block structure for Priority LASSO, which selects features block-wise.
#'   \item Train a Random Forest (\code{ranger}) model on the selected features using \emph{k}-fold cross-validation.
#'   \item Evaluate predictive performance using multiple regression metrics (e.g., R-squared, RMSE, MAE, MAPE).
#'   \item Repeat the above for different \code{k} values and clustering methods, and identify the best-performing configuration.
#' }
#'
#' The primary objective is to compare the performance of different clustering methods
#' when used to define blocks for Priority LASSO, and to select the best-performing
#' clustering method as the final recommended approach.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{<method>}{For each clustering method used, a list containing:
#'     \describe{
#'       \item{aggregates}{A data frame of aggregated performance metrics for each \code{k} tested.}
#'       \item{scores}{A data frame of per-fold cross-validation scores for each \code{k} tested.}
#'     }
#'   }
#'   \item{model_best}{A data frame summarizing the best-performing \code{k} for each clustering method, with the selected features and best R-squared.}
#'   \item{feature_selection}{A character vector of selected features from the best overall model.}
#'   \item{time_mins}{Total computation time in minutes.}
#'   \item{memory_bytes}{Total memory usage in bytes during computation.}
#' }
#'
#' @references
#' Klau S, Jurinovic V, Hornung R, et al (2018) Priority-Lasso: a simple hierarchical approach to the prediction of clinical outcome using multi-omics data. BMC bioinformatics 19(1):322.
#' \url{https://doi.org/10.1186/s12859-018-2344-6}
#'
#' Kaufman L, Rousseeuw PJ (2009) Finding groups in data: an introduction to cluster analysis. John Wiley & Sons. ISBN:9780470317488.
#' \href{https://books.google.co.uk/books?hl=zh-CN&lr=&id=YeFQHiikNo0C&oi=fnd&pg=PR11&dq=Kaufman+L,+Rousseeuw+PJ+(2009)+Finding+groups+in+data:+an+introduction+to+cluster+analysis.+John+Wiley+%26+Sons.&ots=5CscG0KDxz&sig=4DInpEYaSlQcdQm3jlcJ0YaCJ70&redir_esc=y#v=onepage&q=Kaufman%20L%2C%20Rousseeuw%20PJ%20(2009)%20Finding%20groups%20in%20data%3A%20an%20introduction%20to%20cluster%20analysis.%20John%20Wiley%20%26%20Sons.&f=false}{books.google}
#'
#' MacQueen J (1967) Multivariate observations. Proceedings ofthe 5th Berkeley Symposium on Mathematical Statisticsand Probability 1:281-297.
#' \href{https://books.google.co.uk/books?hl=zh-CN&lr=&id=IC4Ku_7dBFUC&oi=fnd&pg=PA281&dq=MacQueen,+J.+(1967).+Some+methods+for+classification+and+analysis+of+multivariate+observations.+In+Proceedings+of+the+Fifth+Berkeley+Symposium+on+Mathematical+Statistics+and+Probability,+1,+281%E2%80%93297.&ots=nQXfKVGdpK&sig=HIm1jb6jFQt2HIdFxTgW-993Q3E&redir_esc=y#v=onepage&q&f=false}{books.google}
#'
#' Reynolds D (2015) Gaussian mixture models. Encyclopedia of biometrics. Springer, Boston 827-832.
#' \url{https://doi.org/10.1007/978-1-4899-7488-4_196}
#'
#' @examples
#' \dontrun{
#' PL_result <- Cluster_priority_lasso(
#' data = WW1768,
#' cluster_methods = c("clara", "kmeans", "gmm"),
#' k_values = 2:5,
#' clara.samples = 5,
#' seed = 1,
#' n_folds = 5,
#' output_file = "Cluster_priority_lasso.xlsx"
#' )
#' }
#'
#' @seealso \code{\link{boruta_regression},\link{auto_abess},\link{glmnet_regression},\link{mmpc_regression}}}
#'
#' @export
Cluster_priority_lasso <- function(data,
                                   cluster_methods = c("clara", "kmeans", "gmm"),
                                   k_values = NULL, clara.samples = 5, seed = 1, n_folds = 10,
                                   output_file = "Cluster_priority_lasso.xlsx") {
  suppressPackageStartupMessages({
    require(mclust)})

  `%>%` <- magrittr::`%>%`
  start_time <- Sys.time(); gc(); start_mem <- pryr::mem_used()

  all_results <- list()
  metrics <- list(
    rsq = mlr3::msr("regr.rsq"),
    mse = mlr3::msr("regr.mse"),
    time_both = mlr3::msr("time_both"),
    rmse = mlr3::msr("regr.rmse"),
    bias = mlr3::msr("regr.bias"),
    ktau = mlr3::msr("regr.ktau"),
    mae = mlr3::msr("regr.mae"),
    mape = mlr3::msr("regr.mape"),
    maxae = mlr3::msr("regr.maxae"),
    medae = mlr3::msr("regr.medae"),
    medse = mlr3::msr("regr.medse"),
    msle = mlr3::msr("regr.msle"),
    pbias = mlr3::msr("regr.pbias"),
    rae = mlr3::msr("regr.rae"),
    rmsle = mlr3::msr("regr.rmsle"),
    rrse = mlr3::msr("regr.rrse"),
    rse = mlr3::msr("regr.rse"),
    sae = mlr3::msr("regr.sae"),
    smape = mlr3::msr("regr.smape"),
    srho = mlr3::msr("regr.srho"),
    sse = mlr3::msr("regr.sse")
  )

  target <- names(data)[1]
  feature_names <- setdiff(names(data), target)
  features_matrix <- t(data[, feature_names, drop = FALSE])
  rownames(features_matrix) <- feature_names
  n_features <- nrow(features_matrix)

  for (cluster_method in cluster_methods) {
    method_results <- list(
      aggregates = data.frame(),
      scores = data.frame()
    )

    current_k_values <- if (is.null(k_values)) {
      switch(cluster_method,
             "gmm" = 2:min(10, n_features-1),
             1:min(10, n_features-1))
    } else {
      k_values[k_values >= ifelse(cluster_method == "gmm", 2, 1) &
                 k_values <= n_features]
    }

    results_summary <- data.frame()
    fold_results_list <- list()

    for (k in current_k_values) {
      message(paste("\nRunning", cluster_method, "with k =", k))
      cluster_result <- tryCatch({
        switch(cluster_method,
               "clara" = cluster::clara(features_matrix, k = k, samples = clara.samples),
               "kmeans" = stats::kmeans(features_matrix, centers = k),
               "gmm" = mclust::Mclust(features_matrix, G = k))
      }, error = function(e) {
        message(paste(cluster_method, "k =", k, "failed:", e$message))
        return(NULL)
      })

      if (is.null(cluster_result)) next

      clusters <- switch(cluster_method,
                         "clara" = cluster_result$clustering,
                         "kmeans" = cluster_result$cluster,
                         "gmm" = cluster_result$classification)

      blocks_param <- split(names(clusters), clusters) %>%
        lapply(function(group) {
          which(feature_names %in% group)
        })

      task <- mlr3::as_task_regr(data, target = target)
      flt <- mlr3filters::flt("selected_features",
                 learner = mlr3::lrn("regr.priority_lasso", lambda.type = "lambda.min"),
                 blocks = blocks_param)

      tryCatch({
        flt$calculate(task)
        score_vec <- flt$scores
        selected_features <- names(score_vec)[which(score_vec == 1)]

        if (length(selected_features) == 0) {
          message(paste(cluster_method, "k =", k, "no features selected"))
          next
        }
      }, error = function(e) {
        message(paste(cluster_method, "k =", k, "feature selection failed:", e$message))
        return(NULL)
      })

      PL <- data[, c(target, selected_features), drop = FALSE]
      PL_task <- mlr3::as_task_regr(PL, target = target)

      set.seed(seed)
      learner <- mlr3::lrn("regr.ranger")

      rr <- mlr3::resample(PL_task, learner, mlr3::rsmp("cv", folds = n_folds), store_models = FALSE)
      get_metric <- function(m) {rr$aggregate(m)[[1]]}
      metric_values <- base::sapply(metrics, get_metric)

      results_summary <- base::rbind(results_summary, data.frame(
        algorithm = cluster_method,
        k = k,
        rsq = metric_values["rsq"],
        mse = metric_values["mse"],
        time_both = metric_values["time_both"],
        rmse = metric_values["rmse"],
        bias = metric_values["bias"],
        ktau = metric_values["ktau"],
        mae = metric_values["mae"],
        mape = metric_values["mape"],
        maxae = metric_values["maxae"],
        medae = metric_values["medae"],
        medse = metric_values["medse"],
        msle = metric_values["msle"],
        pbias = metric_values["pbias"],
        rae = metric_values["rae"],
        rmsle = metric_values["rmsle"],
        rrse = metric_values["rrse"],
        rse = metric_values["rse"],
        sae = metric_values["sae"],
        smape = metric_values["smape"],
        srho = metric_values["srho"],
        sse = metric_values["sse"],
        num_selected_features = length(selected_features),
        selected_features = paste(selected_features, collapse = ",")
      ))

      fold_results <- tryCatch({
        res <- rr$score(metrics)
        res$algorithm <- cluster_method
        res$k <- k
        res <- res[, !names(res) %in% c("task", "learner", "resampling", "prediction_test")] %>%
          dplyr::select(c("iteration","algorithm","k","task_id",'task_id'),dplyr::everything())
        res
      }, error = function(e) NULL)

      if (!is.null(fold_results)) {
        fold_results_list[[paste(cluster_method, k)]] <- fold_results
      }
    }

    method_results$aggregates <- results_summary
    method_results$scores <- do.call(rbind, fold_results_list)

    all_results[[cluster_method]] <- method_results
  }

  end_time <- Sys.time(); gc(); end_mem <- pryr::mem_used()
  total_time <- round(as.numeric(difftime(end_time, start_time, units = "mins")), 3)
  total_bytes <- end_mem - start_mem

  if (length(all_results) > 0) {
    best_models <- purrr::imap_dfr(all_results, ~ {
      summary_df <- .x$aggregates
      summary_df %>%
        dplyr::filter(rsq == max(rsq, na.rm = TRUE)) %>%
        dplyr::slice(1) %>%
        dplyr::transmute(
          algorithm = .y,
          best_rsq = rsq,
          Selected_Features = selected_features,
          Optimal_k = k
        )
    })

    if (nrow(best_models) > 0) {
      all_results[["model_best"]] <- best_models %>% { `rownames<-`(., NULL) }
    }
  }

  final_output <- list()
  for (method in setdiff(names(all_results), "model_best")) {
    if (nrow(all_results[[method]]$aggregates) > 0) {
      final_output[[paste(method, "aggregates")]] <- all_results[[method]]$aggregates
    }
    if (nrow(all_results[[method]]$scores) > 0) {
      final_output[[paste(method, "scores")]] <- all_results[[method]]$scores
    }
  }

  best_summary <- all_results[["model_best"]]
  best_row <- best_summary[which.max(best_summary$best_rsq), ]
  all_results[["feature_selection"]] <- unlist(strsplit(best_row$Selected_Features, ","))
  all_results[["time_mins"]] <- total_time
  all_results[["memory_bytes"]] <- total_bytes

  final_output[["time_memory"]] <- tibble::tibble(
    time_mins = total_time,memory_bytes = total_bytes)

  writexl::write_xlsx(final_output, path = output_file)
  message(paste("Results saved to:", output_file))

  invisible(all_results)
}
