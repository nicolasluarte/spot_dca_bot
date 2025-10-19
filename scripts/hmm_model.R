# libs ----
pacman::p_load(
    tidyquant,
    ggplot2,
    tidyverse,
    GA,
    caret,
    cluster
)

# data ----
cat("fetching sp500 data...\n")
sp500_data_raw <- tq_get(
    "^GSPC",
    from = "2020-01-01",
    to = Sys.Date()
)
prices <- log(sp500_data_raw$adjusted)
cat("data loaded...\n")
print(head(sp500_data_raw))

# cluster eval ----
evaluate_clustering_cv <- function(window_matrix, num_clusters, k_folds = 3, seed = 123) {
    if (num_clusters <= 1) {
        return(NA)
    }

    set.seed(seed)
    folds <- createFolds(1:nrow(window_matrix), k = k_folds)
    silhouette_scores <- numeric(k_folds)

    for (i in 1:k_folds) {
        train_indices <- unlist(folds[-i])
        test_indices <- folds[[i]]

        train_data <- window_matrix[train_indices, ]
        test_data <- window_matrix[test_indices, ]

        # guard against too many suggested clusters
        num_distinct_points <- nrow(unique(train_data))
        if (num_distinct_points < num_clusters) {
            silhouette_scores[i] <- NA
            next
        }

        kmeans_train <- kmeans(train_data, centers = num_clusters, nstart = 10)

        assign_cluster <- function(point, centroids) {
            distances <- apply(centroids, 1, function(c) sum((point - c)^2))
            return(which.min(distances))
        }

        test_assignments <- apply(test_data, 1, assign_cluster, centroids = kmeans_train$centers)

        if (length(unique(test_assignments)) > 1 && nrow(test_data) > 1) {
            sil_result <- silhouette(test_assignments, dist(test_data))
            silhouette_scores[i] <- mean(sil_result[, "sil_width"])
        } else {
            silhouette_scores[i] <- NA
        }
    }

    return(mean(silhouette_scores, na.rm = TRUE))
}

# fitness function ----
fitness_function <- function(params) {
    d <- round(params[1])
    k <- round(params[2])

    # basic checks for invalid parameters
    if (k <= 1 || d >= length(prices) / 2 || d < 3) {
        return(-Inf)
    }

    # create the window matrix based on 'd'
    num_windows <- length(prices) - d + 1
    window_list <- lapply(1:num_windows, function(i) prices[i:(i + d - 1)])

    # help for normalization
    normalize_window <- function(window) {
        if (sd(window) == 0) {
            return(rep(0, length(window)))
        }
        (window - mean(window)) / sd(window)
    }

    normalized_windows <- lapply(window_list, normalize_window)
    window_matrix <- do.call(rbind, normalized_windows)

    # eval clustering quality
    score <- evaluate_clustering_cv(
        window_matrix,
        num_clusters = k,
        k_folds = 3
    )

    # deal with NA
    if (is.na(score)) {
        return(-Inf)
    }

    return(score)
}

# set GA ----
set.seed(42)
message("start GA...")

lower_bounds <- c(window_size = 3, num_clusters = 2)
upper_bounds <- c(window_size = 7, num_clusters = 50)

ga_results <- ga(
    type = "real-valued",
    fitness = fitness_function,
    lower = lower_bounds,
    upper = upper_bounds,
    popSize = 100,
    maxiter = 500,
    pmutation = 0.2,
    pcrossover = 0.7,
    elitism = 2,
    monitor = TRUE
)

message("\nGA seach complete")
print(summary(ga_results))

# optimal params
optimal_params <- round(ga_results@solution[1, ])
names(optimal_params) <- c("window_size", "num_clusters")
cat("\noptimal parameters found:\n")
print(optimal_params)

# re-create window matrix
optimal_d <- optimal_params[1]
optimal_k <- optimal_params[2]
message("recreating the window matrix with optimal window size...")
num_windows <- length(prices) - optimal_d + 1
window_list <- lapply(1:num_windows, function(i) prices[i:(i + optimal_d - 1)])

# normalization
normalize_window <- function(window) {
    if (sd(window) == 0) {
        return(rep(0, length(window)))
    }
    (window - mean(window)) / sd(window)
}
normalized_windows <- lapply(window_list, normalize_window)
final_window_matrix <- do.call(rbind, normalized_windows)

# train model in full dataset
message("\ntraining the final model of full dataset...")
set.seed(42)
final_model_results <- kmeans(
    final_window_matrix,
    centers = optimal_k,
    nstart = 50
)
message("final model training complete")

final_state_sequence <- final_model_results$cluster
print(head(final_state_sequence))

window_matrix_tibble <- as_tibble(set_names(
    cbind(as.data.frame(final_window_matrix), as.data.frame(final_state_sequence)),
    c(c(1:optimal_d), "cluster")
)) %>%
    mutate(id = row_number()) %>%
    pivot_longer(cols = -c(id, cluster), names_to = "T", values_to = "value") %>%
    mutate(
        id = as.factor(id),
        T = as.numeric(T)
    )
window_matrix_tibble

# plot ----
p1 <- window_matrix_tibble %>%
    ggplot(aes(
        T, value,
        group = id
    )) +
    stat_summary(
        fun.data = "mean_se",
        geom = "line",
        aes(group = 1),
        color = "gray70"
    ) +
    stat_summary(
        fun.data = "mean_se",
        geom = "point",
        aes(group = 1),
        size = 3
    ) +
    facet_wrap(~cluster) +
    scale_x_continuous(breaks = c(1:3)) +
    ggpubr::theme_pubr()
p1
