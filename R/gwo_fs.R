#' Grey Wolf Optimization (GWO) and Prior-Guided Grey Wolf Optimization (GWO) for Feature Selection with Random Forest Evaluation
#'
#' This function performs feature selection using either the standard Grey Wolf Optimization (GWO) algorithm
#' or a Prior-Guided Grey Wolf Optimization (PGWO) variant. In both cases, each candidate solution ("wolf")
#' is represented as a binary vector indicating selected features.
#'
#' @param data A data frame where the first column is the target variable and the remaining columns are SNP features.
#' Feature encoding can be 0/1/2 or -1/0/1.
#' @param feature_cols Optional character vector specifying which columns in \code{data} to use as features.
#' Defaults to all columns except the first (target) column.
#' @param control A list of GWO/PGWO control parameters (population size, max iterations, elite count, etc.).
#' For PGWO, include \code{suggestions} = prior knowledge matrix.
#' Defaults: population size = 50, max iterations = 200, no prior suggestions.
#' @param seed Numeric random seed for reproducibility during Random Forest training.
#' @param n_folds Number of folds for cross-validation (default is 10).
#' @param use_parallel Logical, whether to use parallel processing (default is FALSE).
#' @param n_workers Number of parallel workers; defaults to number of detected CPU cores minus one.
#' @param metrics Character vector of performance metrics to evaluate via \code{mlr3} (default includes R², MSE, RMSE, MAE, etc.).
#' @param output_file Character string specifying the CSV file path to save iteration-wise metrics (default: "GWO.csv").
#'
#' @details
#' The Grey Wolf Optimization algorithm is a swarm intelligence metaheuristic inspired by the social hierarchy
#' and hunting behavior of grey wolves. The population is divided into leaders (alpha, beta, delta) and followers (omega).
#' During each iteration, wolves update their positions based on the positions of the top three leaders, balancing
#' exploration of the feature space with exploitation of promising regions.
#'
#' In the Prior-Guided GWO (PGWO) variant, prior knowledge about potentially relevant features is incorporated via
#' the \code{suggestions} argument in \code{control}. These suggestions are binary feature masks that partially or
#' fully initialize the population, enabling the algorithm to start its search near known informative regions.
#' This can improve convergence speed and solution quality when reliable prior information is available.
#'
#' The fitness function for each wolf is based on the cross-validated performance of a Random Forest regression model
#' (via \code{ranger}), evaluated with metrics provided in the \code{metrics} argument. By default, the algorithm
#' maximizes R², but any set of \code{mlr3} regression measures can be used.
#'
#' Parallel processing can be enabled by setting \code{use_parallel = TRUE}, which distributes fitness evaluations
#' across available CPU cores.
#'
#' @return A list containing:
#' \describe{
#'   \item{performance}{Data frame with detailed iteration-wise performance metrics for all wolves in the population.}
#'   \item{rsq}{The best R² score achieved by the optimization process.}
#'   \item{feature_selection}{Character vector of selected feature names from the best-performing wolf.}
#'   \item{time_mins}{Total runtime of the algorithm in minutes.}
#'   \item{memory_bytes}{Difference in memory usage (in bytes) before and after the run.}
#' }
#'
#' @references
#' Mirjalili S, Mirjalili SM, Lewis A (2014) Grey wolf optimizer. Advances in engineering software 69:46-61.
#' \url{https://doi.org/10.1016/j.advengsoft.2013.12.007}
#'
#' Emary E, Zawbaa HM, Hassanien AE (2016) Binary grey wolf optimization approaches for feature selection. Neurocomputing 172:371-381.
#' \url{https://doi.org/10.1016/j.neucom.2015.06.083}
#'
#' @seealso \code{\link{gwo_plot}}
#'
#' @examples
#' \dontrun{
#' # Standard GWO (random initialization)
#' union <- feature_union(FIAs_W_result$results$FIAs_W, abess_result, Boruta_result, PL_result, EN_result, mmpc_result,data = WW1768)
#' feature_union <- union$union_features
#'
#' GWO_result <- gwo_fs(
#'   data = WW1768,
#'   feature_cols = feature_union,
#'   seed = 1,
#'   control = list(numPopulation = 10, maxIter = 20, elite = 1),
#'   output_file = "GWO.csv"
#' )
#'
#' # Prior-Guided GWO (PGWO) with prior knowledge matrix
#' Prior_matrix <- union$Prior_Guided_matrix
#'
#' PGWO_result <- gwo_fs(
#'   data = WW1768,
#'   feature_cols = feature_union,
#'   seed = 1,
#'   control = list(numPopulation = 10, maxIter = 20, elite = 1, suggestions = Prior_matrix),
#'   output_file = "PGWO.csv"
#' )
#' }
#'
#' @export
gwo_fs <- function(
    data,
    feature_cols = NULL,
    control = list(),
    seed = 1,
    n_folds = 10,
    use_parallel = FALSE,
    n_workers = NULL,
    metrics = c("regr.rsq", "regr.mse", "time_both", "regr.rmse", "regr.bias",
                "regr.ktau", "regr.mae", "regr.mape", "regr.maxae", "regr.medae",
                "regr.medse", "regr.msle", "regr.pbias", "regr.rae", "regr.rmsle",
                "regr.rrse", "regr.rse", "regr.sae", "regr.smape", "regr.srho", "regr.sse"),
    output_file = "GWO.csv") {
  `%>%` <- magrittr::`%>%`
  if (use_parallel) {
    if (is.null(n_workers)) {
      n_workers <- parallel::detectCores() - 1
    }
    future::plan(future::multisession, workers = n_workers)
  }

  start_time <- Sys.time(); gc(); start_mem <- pryr::mem_used()

  if (!is.data.frame(data)) stop("Input data must be a data frame")
  target_col <- names(data)[1]
  message("\nAuto-selected target column: ", target_col)
  if (is.null(feature_cols)) {
    feature_cols <- names(data)[-1]
    message("\nAutomatically detected ", length(feature_cols), " feature columns.")
  }
  missing_cols <- setdiff(feature_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("\nMissing columns detected: ", paste(missing_cols, collapse = ", "))
  }

  X <- data[, feature_cols, drop = FALSE]
  performance <- data.frame()

  fitness_function <- function(x) {
    binary_set <- x > 0.5
    select_col <- which(binary_set)
    set.seed(seed)
    if (length(select_col) == 0) return(0)
    Select_X <- X[, select_col, drop = FALSE]
    GWO_tk <- data %>% dplyr::select(1, names(Select_X))
    task <- mlr3::as_task_regr(GWO_tk, target = target_col)
    learner <- mlr3::lrn("regr.ranger")
    resampling <- mlr3::rsmp("cv", folds = n_folds)
    rr <- mlr3::resample(task, learner, resampling)

    metrics_agg <- rr$aggregate(msrs(metrics))
    metrics_vct <- as.numeric(metrics_agg)
    iteration_count <- nrow(performance) + 1
    features_string <- paste(select_col, collapse = ",")

    new_row <- data.frame(
      iteration = iteration_count,
      particle = features_string,
      as.list(setNames(metrics_vct, metrics)),
      stringsAsFactors = FALSE
    )

    performance <<- rbind(performance, new_row)
    R2 <- rr$aggregate(msrs("regr.rsq"))
    message(sprintf("Iteration: %d | R²: %.5f | Features: %d",
                    iteration_count, R2, length(select_col)))
    return(R2)
  }

  default_control <- list(
    optimType = "MAX",
    numVar = ncol(X),
    numPopulation = 50,
    maxIter = 200,
    rangeVar = matrix(rep(c(0, 1), each = ncol(X)), nrow = 2, byrow = TRUE),
    suggestions = NULL,
    elite = NULL
  )

  final_control <- modifyList(default_control, control)

  gwo_elite <- function(FUN, optimType = "MAX", numVar, numPopulation = 50,
                         maxIter = 200, rangeVar, suggestions = NULL, elite = NULL) {
    dimVar <- numVar
    if (is.null(suggestions)) {
      pop <- matrix(runif(numPopulation * dimVar, rangeVar[1,], rangeVar[2,]),
                    nrow = numPopulation, byrow = TRUE)
    } else {
      suggestions <- as.matrix(suggestions)
      if (ncol(suggestions) != dimVar) {
        stop("Column number of 'suggestions' must match numVar")
      }
      num_suggestions <- nrow(suggestions)
      if (num_suggestions >= numPopulation) {
        pop <- suggestions[1:numPopulation, ]
      } else {
        num_random <- numPopulation - num_suggestions
        random_part <- matrix(runif(num_random * dimVar, rangeVar[1,], rangeVar[2,]),
                              nrow = num_random, byrow = TRUE)
        pop <- rbind(suggestions, random_part)
      }
    }

    fitness <- apply(pop, 1, FUN)
    if (optimType == "MIN") fitness <- -fitness

    alpha_pos <- pop[which.max(fitness), ]
    alpha_score <- max(fitness)

    temp <- fitness; temp[which.max(fitness)] <- -Inf
    beta_pos <- pop[which.max(temp), ]
    beta_score <- max(temp)

    temp[which.max(temp)] <- -Inf
    delta_pos <- pop[which.max(temp), ]
    delta_score <- max(temp)

    for (iter in 1:maxIter) {
      a <- 2 - iter * (2 / maxIter)

      if (!is.null(elite) && elite > 0) {
        elite_indices <- order(fitness, decreasing = TRUE)[1:elite]
        elite_pop <- pop[elite_indices, , drop = FALSE]
        elite_fitness <- fitness[elite_indices]
      }

      for (i in 1:numPopulation) {
        for (j in 1:dimVar) {
          r1 <- runif(1); r2 <- runif(1)
          A1 <- 2 * a * r1 - a; C1 <- 2 * r2
          D_alpha <- abs(C1 * alpha_pos[j] - pop[i,j])
          X1 <- alpha_pos[j] - A1 * D_alpha

          r1 <- runif(1); r2 <- runif(1)
          A2 <- 2 * a * r1 - a; C2 <- 2 * r2
          D_beta <- abs(C2 * beta_pos[j] - pop[i,j])
          X2 <- beta_pos[j] - A2 * D_beta

          r1 <- runif(1); r2 <- runif(1)
          A3 <- 2 * a * r1 - a; C3 <- 2 * r2
          D_delta <- abs(C3 * delta_pos[j] - pop[i,j])
          X3 <- delta_pos[j] - A3 * D_delta

          pop[i,j] <- (X1 + X2 + X3) / 3
        }
        pop[i,] <- pmax(pmin(pop[i,], rangeVar[2,]), rangeVar[1,])
      }

      fitness <- apply(pop, 1, FUN)
      if (optimType == "MIN") fitness <- -fitness

      cat(sprintf("Generation %d: Max = %.6f, Median = %.6f, Mean = %.6f\n",
                  iter, max(fitness), median(fitness), mean(fitness)))

      if (!is.null(elite) && elite > 0) {
        worst_indices <- order(fitness)[1:elite]
        pop[worst_indices, ] <- elite_pop
        fitness[worst_indices] <- elite_fitness
      }

      best_index <- which.max(fitness)
      if (fitness[best_index] > alpha_score) {
        alpha_score <- fitness[best_index]
        alpha_pos <- pop[best_index, ]
      }

      temp <- fitness; temp[which.max(fitness)] <- -Inf
      beta_pos <- pop[which.max(temp), ]
      beta_score <- max(temp)

      temp[which.max(temp)] <- -Inf
      delta_pos <- pop[which.max(temp), ]
      delta_score <- max(temp)
    }

    return(list(best_pos = alpha_pos,
                best_score = if (optimType == "MAX") alpha_score else -alpha_score))
  }

  do.call(gwo_elite, c(list(FUN = fitness_function), final_control))
  end_time <- Sys.time(); gc(); end_mem <- pryr::mem_used()
  total_time <- round(as.numeric(difftime(end_time, start_time, units = "mins")), 3)
  total_bytes <- end_mem - start_mem

  names(performance) <- c("iteration","particle",metrics)

  n_iter <- final_control$maxIter
  n_particle <- final_control$numPopulation
  if (nrow(performance) == (n_iter + 1) * n_particle) {
    performance$iteration <- rep(c(0,seq_len(n_iter)), each = n_particle)
  }

  performance <- performance %>%
    dplyr::mutate(dplyr::across(3:ncol(performance), as.numeric))

  data.table::fwrite(performance, file = output_file)
  message("\nSuccessfully saved GWO results to: ", normalizePath(output_file))

  row <- performance[which.max(performance$regr.rsq), ]
  indices <- as.numeric(strsplit(row$particle, ",")[[1]])
  features <- feature_cols[indices]

  if (use_parallel) {
    future::plan(future::sequential)
  }

  return(list(
    performance = performance,
    rsq = row$regr.rsq,
    feature_selection = features,
    time_mins = total_time,
    memory_bytes = total_bytes
  ))
}
