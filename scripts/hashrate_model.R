pacman::p_load(
    tidyverse,
    ggplot2,
    httr,
    jsonlite,
    randomForest
)
setwd(this.path::here())

hashrate_url <- "https://api.blockchain.info/charts/hash-rate?timespan=all&format=json"
btc_url <- "https://api.blockchain.info/charts/market-price?timespan=all&format=json"
trans_url <- "https://api.blockchain.info/charts/n-transactions?timespan=all&format=json"
fees_url <- "https://api.blockchain.info/charts/transaction-fees-usd?timespan=all&format=json"

url_list <- list(
    hashrate_url,
    btc_url,
    trans_url,
    fees_url
)

data_names <- c(
    "hashrate",
    "btc_price",
    "transactions",
    "fees"
)

data <- url_list %>%
    imap(
        ., function(X, idx) {
            response <- GET(X)
            json_content <- content(response, as = "text")
            data_list <- fromJSON(json_content)
            price_df <- tibble(data_list$values) %>%
                mutate(class = data_names[idx])
            return(price_df)
        }
    )

# fred data
dff <- read_csv("DFF.csv") %>%
    mutate(
        observation_date = as.numeric(as.POSIXct(observation_date, tz = "UTC")),
        class = "dff"
    ) %>%
    rename(
        x = observation_date,
        y = DFF
    )
walcl <- read_csv("WALCL.csv") %>%
    mutate(
        observation_date = as.numeric(as.POSIXct(observation_date, tz = "UTC")),
        class = "walcl"
    ) %>%
    rename(
        x = observation_date,
        y = WALCL
    )


l_n <- 1
data_full <- bind_rows(data, dff, walcl) %>%
    pivot_wider(
        names_from = class,
        values_from = y
    ) %>%
    arrange(x) %>%
    mutate(
        price_int = zoo::na.approx(btc_price, rule = 2),
        hashrate_int = zoo::na.approx(hashrate, rule = 2),
        transactions_int = zoo::na.approx(transactions, rule = 2),
        fees_int = zoo::na.approx(fees, rule = 2),
        dff_int = zoo::na.approx(dff, rule = 2),
        walcl_int = zoo::na.approx(walcl, rule = 2),
        price_lag = lag(price_int, n = l_n),
        hashrate_lag = lag(hashrate_int, n = l_n),
        transactions_lag = lag(transactions_int, n = l_n),
        fees_lag = lag(fees_int, n = l_n),
        dff_lag = lag(dff_int, n = l_n),
        walcl_lag = lag(walcl_int, n = l_n),
        log_return = log(price_int) - log(lag(price_int, n = 1)),
        return_bin = if_else(exp(log_return) - 1 > 0, 1, 0)
    ) %>%
    drop_na(log_return) %>%
    filter(log_return < Inf)
data_mdl <- lm(
    data = data_full,
    log_return ~ log(hashrate_lag + 1) *
        log(transactions_lag + 1) *
        log(fees_lag + 1) *
        log(walcl_lag + 1) *
        log(dff_lag + 1)
)
summary(data_mdl)

data_mdl <- glm(
    data = data_full,
    return_bin ~ log(hashrate_lag + 1) *
        log(transactions_lag + 1) *
        log(fees_lag + 1) *
        log(walcl_lag + 1) *
        log(dff_lag + 1),
    family = binomial(link = "logit")
)
summary(data_mdl)

rf_data <- data_full %>%
    mutate(
        date = as_datetime(x),
        return_bin = as.factor(return_bin),
        log_hashrate_lag = log(hashrate_lag + 1),
        log_transactions_lag = log(transactions_lag + 1),
        log_fees_lag = log(fees_lag + 1),
        log_walcl_lag = log(walcl_lag + 1),
        log_dff_lag = log(dff_lag + 1)
    )

rf_mdl <- randomForest(
    return_bin ~ log_hashrate_lag +
        log_transactions_lag +
        log_fees_lag +
        log_walcl_lag +
        log_dff_lag,
    data = rf_data
)

pred_data <- rf_data %>%
    mutate(
        .pred = predict(rf_mdl)
    )

nd <- tail()



probabilities <- predict(data_mdl, type = "response")
predicted_classes <- ifelse(probabilities > 0.5, 1, 0)
conf_matrix <- table(Predicted = predicted_classes, Actual = data_full$return_bin)
conf_matrix
accuracy <- (conf_matrix[1, 1] + conf_matrix[2, 2]) / sum(conf_matrix)
print(paste("Accuracy:", accuracy))
