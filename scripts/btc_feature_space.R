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
    tseries,
    ggthemes,
    mclust
)
pacman::p_unload(MASS)
plan(multisession)

setwd(this.path::here())
palette("Okabe-Ito")

# download btc data ----

btc_data <- crypto2::crypto_list() %>%
    filter(slug == "bitcoin") %>%
    crypto2::crypto_history()

# lagged predictors ----
volatility_data <- btc_data %>%
    mutate(
        bin_profit = if_else(close > open, 1, 0),
        lag1 = lag(x = bin_profit, n = 1),
        lag7 = lag(x = bin_profit, n = 7),
        lag14 = lag(x = bin_profit, n = 14),
        lag21 = lag(x = bin_profit, n = 21),
        lag28 = lag(x = bin_profit, n = 28),
        lag35 = lag(x = bin_profit, n = 35)
    ) %>%
    drop_na()

# logistic model ----
l_mdl <- glm(
    data = volatility_data %>%
        select(bin_profit, contains("lag")),
    bin_profit ~ .,
    family = binomial(link = "logit")
)
summary(l_mdl)

l_emm <- emmeans::emmeans(
    l_mdl,
    specs = ~.,
    type = "response"
)
l_emm
