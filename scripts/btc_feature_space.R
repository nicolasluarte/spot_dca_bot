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

## simple residual model ----

mdl_data <- btc_data %>%
    mutate(
        log_open = log(open),
        log_close = log(close),
        log_high = log(high),
        log_low = log(low),
        t = scales::rescale(timestamp, to = c(0, 1))
    )
mdl_data

mdl_open <- lm(
    data = mdl_data,
    log_open ~ t
)
mdl_close <- lm(
    data = mdl_data,
    log_close ~ t
)
mdl_high <- lm(
    data = mdl_data,
    log_high ~ t
)
mdl_low <- lm(
    data = mdl_data,
    log_low ~ t
)

sigma(mdl_open)
sigma(mdl_close)
sigma(mdl_high)
sigma(mdl_low)
