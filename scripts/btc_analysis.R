# libs ----
pacman::p_load(
    tidyverse,
    crypto2,
    furrr,
    fitdistrplus,
    mgcv,
    patchwork,
    extrafont,
    tidymodels,
    ranger
)

setwd(this.path::here())
palette("Okabe-Ito")

# functions ----

# compute entropy per ID per rel_date per exp_group per exp_phase
binary_entropy <- function(p) {
    if (p == 0 || p == 1) {
        return(0)
    }
    return(-(p * log2(p) + (1 - p) * log2(1 - p)))
}


# download btc data ----

btc_data <- crypto2::crypto_list() %>%
    filter(slug == "bitcoin") %>%
    crypto2::crypto_history()

# make features ----

## coefficient of variation ----

feature_data <- btc_data %>%
    rowwise() %>%
    mutate(
        mean_daily = mean(c_across(c("open", "high", "low", "close"))),
        sd_daily = sd(c_across(c("open", "high", "low", "close")))
    ) %>%
    ungroup() %>%
    mutate(
        cv_daily = sd_daily / mean_daily,
        scaled_ts = scales::rescale(as.numeric(timestamp)),
        detrend_close = scales::rescale(residuals(lm(log(close) ~ scaled_ts)), to = c(0, 1))
    )

# fit parametric distribution to cv

fit_g <- fitdist(feature_data$cv_daily,
    "gamma",
    method = "mle"
)

plot(fit_g)

# daily coefficient of variation parametric fit

parametric_cv <- feature_data %>%
    mutate(
        gamma_fit = pgamma(cv_daily,
            shape = fit_g$estimate[1],
            rate = fit_g$estimate[2]
        )
    )

# models ----

## original data ----

rf_data <- parametric_cv %>%
    arrange(scaled_ts) %>%
    select(
        timestamp,
        detrend_close,
        scaled_ts,
        gamma_fit,
        sd_daily,
        cv_daily,
        open,
        high,
        low,
        close
    )

## training and test data -----
set.seed(420)
data_split <- initial_time_split(
    rf_data,
    prop = 0.8
)
train_raw_data <- training(data_split)
test_raw_data <- testing(data_split)

## feature engineering ----
lag_recipe <- recipe(
    detrend_close ~ .,
    data = train_raw_data
) %>%
    step_lag(
        timestamp,
        gamma_fit,
        sd_daily,
        cv_daily,
        open,
        high,
        low,
        close,
        lag = 1:30,
        prefix = "lag_"
    ) %>%
    step_date(
        timestamp,
        features = c("month", "dow", "doy")
    ) %>%
    step_normalize(all_numeric_predictors()) %>%
    step_naomit(everything())

## random forest specification ----

rf_model_spec <- rand_forest(
    mtry = tune(),
    trees = 500,
    min_n = tune()
) %>%
    set_engine("ranger", importance = "permutation") %>%
    set_mode("regression")

## workflow ----

rf_workflow <- workflow() %>%
    add_recipe(lag_recipe) %>%
    add_model(rf_model_spec)

## time series cv ----
set.seed(69)
time_series_cv_folds <- rolling_origin(
    train_raw_data,
    initial = floor(nrow(train_raw_data) * 0.6),
    asses = floor(nrow(train_raw_data) * 0.15),
    skip = floor(nrow(train_raw_data) * 0.10),
    cumulative = FALSE
)

## hyperparameter tuning grid ----
temp_prepped_recipe <- prep(
    lag_recipe,
    training = train_raw_data
)

num_predictor_approx <- length(
    temp_prepped_recipe$term_info %>%
        filter(role == "predictor") %>%
        pull(variable)
)

mtry_upper_bound <- max(2, num_predictor_approx)
mtry_lower_bound <- 2

rf_param_grid <- grid_latin_hypercube(
    mtry(range = c(mtry_lower_bound, mtry_upper_bound)),
    min_n(range = c(2, 25)),
    size = 100
)

## tunning process ----
set.seed(43)
control_tuning <- control_grid(
    save_pred = TRUE,
    verbose = TRUE,
    parallel_over = "everything"
)

rf_tune_results <- tune_grid(
    rf_workflow,
    resamples = time_series_cv_folds,
    grid = rf_param_grid,
    metrics = metric_set(rmse, mae, rsq),
    control = control_tuning
)

## analyze tunning results ----

show_best(
    rf_tune_results,
    metric = "rmse",
    n = 5
)

autoplot(rf_tune_results)

best_rf_params <- select_best(
    rf_tune_results,
    metric = "rmse"
)
best_rf_params

## finalize workflow ----

final_rf_workflow <- finalize_workflow(
    rf_workflow,
    best_rf_params
)

final_rf_model_fit <- fit(
    final_rf_workflow,
    data = train_raw_data
)

## evaluate final model ----

processed_test_data_for_eval <- prep(
    lag_recipe,
    training = train_raw_data
) %>%
    bake(new_data = test_raw_data)

augment_test_data <- augment(
    final_rf_model_fit,
    new_data = test_raw_data
)

final_metrics <- yardstick::metrics(
    augment_test_data,
    truth = detrend_close,
    estimate = .pred
)
final_metrics

augment_test_data %>%
    filter(timestamp > "2024-01-01") %>%
    ggplot(aes(
        x = timestamp
    )) +
    geom_line(aes(
        y = detrend_close,
        color = "actual"
    )) +
    geom_line(aes(
        y = .pred,
        color = "predicted"
    ), linetype = "dashed")
