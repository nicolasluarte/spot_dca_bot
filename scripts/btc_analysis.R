# libs ----
pacman::p_load(
    tidyverse,
    crypto2,
    furrr,
    fitdistrplus,
    patchwork,
    extrafont,
    tidymodels,
    ranger,
    future,
    TTR,
    janitor,
    vip,
    tseries
)
pacman::p_unload(MASS)
plan(multisession)

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

RSI_trend <- TTR::RSI(price = parametric_cv$close, n = 14)
ATR_trend <- TTR::ATR(HLC = parametric_cv[, c("high", "low", "close")], n = 14)
BB_trend <- TTR::BBands(HLC = parametric_cv[, c("high", "low", "close")], n = 14)
ROC_trend <- TTR::ROC(x = parametric_cv$close, n = 9)
STOCH_trend <- TTR::stoch(HLC = parametric_cv[, c("high", "low", "close")])
ADX_trend <- TTR::ADX(HLC = parametric_cv[, c("high", "low", "close")], n = 14)

rf_data <- parametric_cv %>%
    arrange(scaled_ts) %>%
    select(
        timestamp,
        detrend_close,
        gamma_fit,
        sd_daily,
        cv_daily,
        open,
        high,
        low,
        close
    ) %>%
    mutate(
        RSI_trend = as_tibble(RSI_trend)$value,
        ATR_trend = as_tibble(ATR_trend)$atr,
        BB_trend = as_tibble(BB_trend)$pctB,
        ROC_trend = as_tibble(ROC_trend)$value,
        STOCH_fastK = as_tibble(STOCH_trend)$fastK,
        STOCH_fastD = as_tibble(STOCH_trend)$fastD,
        STOCH_slowD = as_tibble(STOCH_trend)$slowD,
        ADX_adx = as_tibble(ADX_trend)$ADX,
        ADX_pos = as_tibble(ADX_trend)$DIp,
        ADX_neg = as_tibble(ADX_trend)$DIn
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
        RSI_trend,
        ATR_trend,
        BB_trend,
        ROC_trend,
        STOCH_fastK,
        STOCH_fastD,
        STOCH_slowD,
        ADX_adx,
        ADX_pos,
        ADX_neg,
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
    trees = 500
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
    initial = floor(nrow(train_raw_data) * 0.3),
    asses = floor(nrow(train_raw_data) * 0.05),
    skip = floor(nrow(train_raw_data) * 0.05),
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

rf_param_grid <- grid_space_filling(
    mtry(range = c(mtry_lower_bound, mtry_upper_bound)),
    size = 25
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

### variable importance ----
importance_plot <- vip(
    final_rf_model_fit,
    num_features = 50
)
importance_plot


residual_data <- augment_test_data %>%
    mutate(
        y = scales::rescale(detrend_close, to = c(0, 1)),
        y_hat = scales::rescale(.pred, to = c(0, 1)),
        .resid = y - y_hat
    ) %>%
    select(
        timestamp,
        .resid
    ) %>%
    mutate(
        x_t = lag(.resid, n = 1),
        delta_x = (.resid - x_t)
    )

adf.test(residual_data$.resid)

plot(density(residual_data$.resid))

residual_data %>%
    drop_na() %>%
    ggplot(aes(
        x_t, delta_x
    )) +
    geom_point() +
    geom_smooth(method = "lm") +
    facet_wrap(~ x_t > 0, scales = "free_x")

augment_test_data %>%
    ggplot(aes(
        x = timestamp
    )) +
    geom_line(aes(
        y = scales::rescale(detrend_close, to = c(0, 1)),
        color = "actual"
    )) +
    geom_line(aes(
        y = scales::rescale(.pred, to = c(0, 1)),
        color = "predicted"
    )) +
    scale_color_manual(
        values = c("black", "orange")
    ) +
    ggpubr::theme_pubr() +
    labs(
        color = ""
    ) +
    ylab("Detrended rescaled residuals")
