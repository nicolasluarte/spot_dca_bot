pacman::p_load(
    tidyverse,
    ggplot2,
    furrr,
    nloptr,
    crypto2,
    splines,
    Rcpp
)

# Terminal args ----

# note that this is a vector
args <- commandArgs(trailingOnly=TRUE)

# print coin
print(paste("Selected coin:", args[1], sep = " "))

# https://stackoverflow.com/questions/47044068/get-the-path-of-current-script
# get path of source file
getCurrentFileLocation <-  function()
{
    this_file <- commandArgs() %>% 
        tibble::enframe(name = NULL) %>%
        tidyr::separate(col=value, into=c("key", "value"), sep="=", fill='right') %>%
        dplyr::filter(key == "--file") %>%
        dplyr::pull(value)
    if (length(this_file)==0)
    {
        this_file <- rstudioapi::getSourceEditorContext()$path
    }
    return(dirname(this_file))
}
# sets path on source file location
script_path <- getCurrentFileLocation()
setwd(script_path)

# get crypto data ----

# create dir to store all crypto data
dir.create("data")
dir.create("optimization_results")

# get list of all relevant crypto coins
if (!file.exists("data/crypto_ls.rds")){
    print("Creating crypto list...")
    crypto_ls <-
        crypto2::crypto_list() %>% 
        filter(name %in% c("Bitcoin", "NEAR Protocol", "PAX Gold"))
    write_rds(crypto_ls, "data/crypto_ls.rds")
} else{
    print("Crypto list already exists...")
    crypto_ls <- read_rds("data/crypto_ls.rds")
}

# download list of selected crypto coins
if(!file.exists("data/crypto_data.rds")){
    print("Downloading crypto data...")
    crypto_data <-
        crypto2::crypto_history(
            read_rds("data/crypto_ls.rds"),
            convert = "USD"
        )
    write_rds(crypto_data, "data/crypto_data.rds")
} else{
    print("Crypto data already downloaded...")
    crypto_data <- read_rds("data/crypto_data.rds")
}

# Simulation definitions ----
coin <- args[1]
n_days <- as.numeric(args[2])
n_boots <- as.numeric(args[3])
n_simulations <- as.numeric(args[4])
CAPITAL <- 10000


# functions ----

extract_nbars <-
    function(x, len){
        sample_day <- sample(head(row_number(x), -len), size = 1, replace = TRUE)
        filtered_data <-
            x %>% 
            mutate(idx = row_number()) %>% 
            filter(
                idx >= sample_day,
                idx <= sample_day + (len-1)
            ) 
        return(filtered_data)
    }

do_linear_reg <-
    function(x, p){
        x <- x %>% 
            mutate(
                close = close/(close[1]),
            )
        mdl <-
            lm(
                data = x,
                close ~ poly(as.numeric(as.factor(timestamp)),
                             degree = 3,
                             raw = TRUE)
            )
        tidy_lm <-
            mdl %>% 
            broom::tidy()
        stats <-
            broom::glance(mdl)
        if (p == "plot"){
            p1 <-
                x %>% 
                ggplot(aes(
                    as.numeric(as.factor(timestamp)), (close)
                )) +
                geom_line() +
                geom_point() +
                geom_smooth(method = "lm", se = FALSE)
                theme_classic()
            return(list(mdl = mdl, mdl_tidy = tidy_lm, stats=stats, plot = p1))
        }
        else{
            return(list(mdl = mdl, mdl_tidy=tidy_lm, stats=stats))
        }
    }

bootstrap_stats <-
    function(nboot, len, data){
        bootstraps <-
            1:nboot %>% 
            map_dfr(
                ., function(x){
                    d <-extract_nbars(data, len) %>% 
                        do_linear_reg(., "")
                    intercept <- d$mdl_tidy$estimate[1]
                    poly_1 <- d$mdl_tidy$estimate[2]
                    poly_2 <- d$mdl_tidy$estimate[3]
                    poly_3 <- d$mdl_tidy$estimate[4]
                    mdl_sigma <- d$stats$sigma[1]
                    return(tibble(
                        intercept = intercept,
                        poly_1 = poly_1,
                        poly_2 = poly_2,
                        poly_3 = poly_3,
                        mdl_sigma = mdl_sigma
                    ))
                }
            )
        return(bootstraps)
    }

run_sim <- function(jds, len, mod){
    i <- 1:(len)
    slope <- (
        (jds$poly_1 * (i)) + 
            (jds$poly_2 * (i)^2) +
            (jds$poly_3 * (i)^3)
        ) + jds$intercept
    rand_noise <-
        rnorm(
            mean = slope,
            sd = jds$mdl_sigma,
            n = len
        )
    # normalize to avoid negative values, and only get the 
    # shape of the price,, I add a constant to avoid price 0
    # that would lead to an infinite buy
    # we take the first one out in case is 0
    d <- (slope + rand_noise)
    return(tibble(
        x = 1:(len * mod),
        y = scales::rescale(d, to=c(1,2))
    ))
}

# bootstraps ----

# some sample have HUGE sigma, so I just filtered them out
coin_data <- crypto_data %>%
    filter(symbol == coin)
print(paste("Dimension of data", dim(coin_data), sep = " "))
print("Bootstrapping samples...")
b <- bootstrap_stats(n_boots, n_days, coin_data) %>% 
    drop_na() %>% 
    filter(mdl_sigma < 1)

# multivariate distribution sampling ----

covariance_matrix <- cov(b)
mean_matrix <- c(
    median(b$intercept),
    median(b$poly_1),
    median(b$poly_2),
    median(b$poly_3),
    median(b$mdl_sigma)
)

sample_distribution <- MASS::mvrnorm(
    n = n_boots,
    mu = mean_matrix,
    Sigma = covariance_matrix) %>% 
    as_tibble() %>% 
    filter(mdl_sigma >= 0, intercept >= 0)
sample_distribution

# simulations ----

# this simulation is a tick every 10 minutes over 30 days
print("Simulating runs...")
runs <-
    sample_distribution %>% 
    mutate(r = row_number()) %>% 
    group_by(r) %>% 
    group_split() %>% 
    map_dfr(
        ., function(X){
            # 90 days in seconds
            sim <- run_sim(X, n_days, 1) %>% 
                mutate(
                    iteration = X$r,
                    intercept = X$intercept,
                    poly_1 = X$poly_1,
                    poly_2 = X$poly_2,
                    poly_3 = X$poly_3,
                    mdl_sigma = X$mdl_sigma
                    )
            return(sim)
        },
        .progress = TRUE
    )

# DCA loop

# cost function ----

cost_function <- function(
        price_series, # a given simulation
        base_order_size,
        DCA_order_size,
        num_of_dca, # this is the percent of max possible
        dca_multiplier, # multiply each successive dca order
        price_deviation_target, # when to make a dca order
        price_deviation_multiplier,
        take_profit_target, # when to take profit
        days,
        len
){
    # initialize parameters
    days <- days
    len <- days
    price_series <- price_series
    capital <- CAPITAL
    num_of_dca <- num_of_dca
    dca_multiplier <- dca_multiplier
    price_deviation_target <- price_deviation_target
    take_profit_target <- take_profit_target
    dca_counter <- 0
    dca_orders_max <- 1
    average_buy_price <- 0
    price_deviation <- 0
    drawdown <- 0
    drawdown_out <- 0
    
    # simulation parameters
    coin_bag <- c(0) # stores bought coins
    buy_prices <- c(0) # stores coin buy prices
    outflow <- c(0) # stores money spent
    realized_profit <- c(0) # stores realized profit
    DCA_bool <- FALSE # bool to store event
    TP_bool <- FALSE # bool to store event
    
    # create tibble to store simulation variables
    headers <- tibble(
        price = 0,
        TP_bool = FALSE,
        DCA_bool = FALSE,
        dca_orders_max = 0,
        coin_bag = 0,
        buy_prices = 0,
        outflow = 0,
        average_buy_price = 0,
        realized_profit = 0,
        price_deviation = 0,
        dca_counter = 0,
        price_deviation_target = price_deviation_target,
        step = 0
    )
    
    # main simulation loop
    for (i in 1:len){
        # observe current price
        price <- price_series$y[i]
    	drawdown <- append(drawdown, if_else(average_buy_price != 0,
				(price - average_buy_price) / price,
				0))
        
        # update headers
        headers <- 
            headers %>% 
            add_row(
            price = price,
            TP_bool = TP_bool,
            DCA_bool = DCA_bool,
            dca_orders_max = dca_orders_max,
            coin_bag = tail(coin_bag, n=1),
            buy_prices = tail(buy_prices, n=1),
            outflow = tail(outflow, n=1),
            average_buy_price = tail(average_buy_price, n=1),
            realized_profit = tail(realized_profit, n=1),
            price_deviation = price_deviation,
            dca_counter = dca_counter,
            price_deviation_target = price_deviation_target,
            step = i
        )
        
        # reset events
        DCA_bool <- FALSE
        TP_bool <- FALSE
        
        # this is the first buy event
        if (i == 1){
            # check maximum number of DCA orders
            dca_sizes <- (DCA_order_size * (dca_multiplier^(0:100)))
            dca_depth <- cumsum(DCA_order_size * (dca_multiplier^(0:100)))
            available_dca_capital <- capital - base_order_size
            max_depth <- max(which(available_dca_capital >= dca_depth))
            dca_orders_max <- trunc(max_depth * num_of_dca)
            coin_bag <- append(coin_bag, base_order_size/price)
            buy_prices <- append(buy_prices, price)
            outflow <- append(outflow, base_order_size)
            average_buy_price <- sum(outflow)/sum(coin_bag)
        }
        
        # compute actual price deviation from previous order
        if (price == 0){
            price_deviation <- 0
            
        }
        else{
            price_deviation <- (price - tail(buy_prices, n=1)) / price
        }
        
        # make dca order?
        if (price_deviation <= price_deviation_target &
            dca_counter <= dca_orders_max){
            dca_counter <- dca_counter + 1
            price_deviation_target <- price_deviation_target *
                (price_deviation_multiplier^dca_counter)
            DCA_bool <- TRUE
            coin_bag <- append(coin_bag, dca_sizes[dca_counter]/price)
            buy_prices <- append(buy_prices, price)
            outflow <- append(outflow, dca_sizes[dca_counter])
            average_buy_price <- sum(outflow)/sum(coin_bag)
        }
        # take profit
        else if ( (price - average_buy_price)/price >= take_profit_target &
                 i > 1) {
            dca_counter <- 0
            TP_bool <- TRUE
	    drawdown_out <- append(drawdown_out, min(drawdown))
	    drawdown <- 0
            realized_profit <- append(realized_profit,
                                      (sum(coin_bag)*price) - sum(outflow)
            )
            # reset all values
            coin_bag <- c()
            buy_price <- c()
            outflow <- c()
            average_buy_price <- 0
            
            # buy ASAP
            coin_bag <- append(coin_bag, base_order_size/price)
            buy_prices <- append(buy_prices, price)
            outflow <- append(outflow, base_order_size)
            average_buy_price <- sum(outflow)/sum(coin_bag)
        }
        # sell all at the end of simulation
        # testing NOT selling all at the end of the simulation
        if (i == len){
            TP_bool <- FALSE
    	    drawdown_out <- append(drawdown_out, min(drawdown))
            # update tibble for the last time
            headers <- 
                headers %>% 
                add_row(
                price = price,
                TP_bool = TP_bool,
                DCA_bool = DCA_bool,
                dca_orders_max = dca_orders_max,
                coin_bag = tail(coin_bag, n=1),
                buy_prices = tail(buy_prices, n=1),
                outflow = tail(outflow, n=1),
                average_buy_price = tail(average_buy_price, n=1),
                realized_profit = tail(realized_profit, n=1),
                price_deviation = price_deviation,
                dca_counter = dca_counter,
                price_deviation_target = price_deviation_target,
                step = i
            )
            # compute end of simulation results
            profit_out <- headers %>% 
                filter(TP_bool == TRUE) %>% 
                pull(realized_profit) %>% 
                {sum(.)}
            number_of_tp_actions <-
                headers %>% 
                filter(TP_bool == TRUE) %>% 
                pull(TP_bool) %>% 
                {length(.)}
            cost_tibble <- headers %>% 
                filter(price != 0, average_buy_price != 0)
            price <- cost_tibble %>% pull(price)
            average_buy_price <- cost_tibble %>% pull(average_buy_price)
            # price and buy price deviations
            # consider this as the sum of squared differences
            # that is, deviation from the price, but we only care 
            # for deviation above price
            dev <- (average_buy_price - price)
            # I take only positives because this are the one we want to reduce
            # when average buy price > price
            # taking al l leads to unwanted behavior
            # take all deviations from price a square them
            dev_positive <- dev[dev>0]^2
            cost <- sum(dev_positive)
        }
        # if (i %% 100 == 0){
        #     print(paste0("progress: ", round((i/len)*100,2), "%" ))
        #     }
    }
    # the lower the better
    return(cost)
}


cost_function2 <- function(
        price_series, # a given simulation
        base_order_size,
        DCA_order_size,
        num_of_dca, # this is the percent of max possible
        dca_multiplier, # multiply each successive dca order
        price_deviation_target, # when to make a dca order
        price_deviation_multiplier,
        take_profit_target, # when to take profit
        days,
        len
){
    # initialize parameters
    days <- days
    len <- 1 * days # to simulate every 10 minutes
    price_series <- price_series
    capital <- CAPITAL
    num_of_dca <- num_of_dca
    dca_multiplier <- dca_multiplier
    price_deviation_target <- price_deviation_target
    take_profit_target <- take_profit_target
    dca_counter <- 0
    dca_orders_max <- 1
    average_buy_price <- 0
    price_deviation <- 0
    drawdown <- 0
    drawdown_out <- 0
    
    # simulation parameters
    coin_bag <- c(0) # stores bought coins
    buy_prices <- c(0) # stores coin buy prices
    outflow <- c(0) # stores money spent
    realized_profit <- c(0) # stores realized profit
    DCA_bool <- FALSE # bool to store event
    TP_bool <- FALSE # bool to store event
    
    # create tibble to store simulation variables
    headers <- tibble(
        price = 0,
        TP_bool = FALSE,
        DCA_bool = FALSE,
        dca_orders_max = 0,
        coin_bag = 0,
        buy_prices = 0,
        outflow = 0,
        average_buy_price = 0,
        realized_profit = 0,
        price_deviation = 0,
        dca_counter = 0,
        price_deviation_target = price_deviation_target,
        step = 0
    )
    
    # main simulation loop
    for (i in 1:len){
        # observe current price
        price <- price_series$y[i]
    	drawdown <- append(drawdown, if_else(average_buy_price != 0,
				(price - average_buy_price) / price,
				0))
        
        # update headers
        headers <- 
            headers %>% 
            add_row(
            price = price,
            TP_bool = TP_bool,
            DCA_bool = DCA_bool,
            dca_orders_max = dca_orders_max,
            coin_bag = tail(coin_bag, n=1),
            buy_prices = tail(buy_prices, n=1),
            outflow = tail(outflow, n=1),
            average_buy_price = tail(average_buy_price, n=1),
            realized_profit = tail(realized_profit, n=1),
            price_deviation = price_deviation,
            dca_counter = dca_counter,
            price_deviation_target = price_deviation_target,
            step = i
        )
        
        # reset events
        DCA_bool <- FALSE
        TP_bool <- FALSE
        
        # this is the first buy event
        if (i == 1){
            # check maximum number of DCA orders
            dca_sizes <- (DCA_order_size * (dca_multiplier^(0:100)))
            dca_depth <- cumsum(DCA_order_size * (dca_multiplier^(0:100)))
            available_dca_capital <- capital - base_order_size
            max_depth <- max(which(available_dca_capital >= dca_depth))
            dca_orders_max <- trunc(max_depth * num_of_dca)
            coin_bag <- append(coin_bag, base_order_size/price)
            buy_prices <- append(buy_prices, price)
            outflow <- append(outflow, base_order_size)
            average_buy_price <- sum(outflow)/sum(coin_bag)
        }
        
        # compute actual price deviation from previous order
        if (price == 0){
            price_deviation <- 0
            
        }
        else{
            price_deviation <- (price - tail(buy_prices, n=1)) / price
        }
        
        # make dca order?
        if (price_deviation <= price_deviation_target &
            dca_counter <= dca_orders_max){
            dca_counter <- dca_counter + 1
            price_deviation_target <- price_deviation_target *
                (price_deviation_multiplier^dca_counter)
            DCA_bool <- TRUE
            coin_bag <- append(coin_bag, dca_sizes[dca_counter]/price)
            buy_prices <- append(buy_prices, price)
            outflow <- append(outflow, dca_sizes[dca_counter])
            average_buy_price <- sum(outflow)/sum(coin_bag)
        }
        # take profit
        else if ( (price - average_buy_price)/price >= take_profit_target &
                 i > 1) {
            dca_counter <- 0
            TP_bool <- TRUE
	    drawdown_out <- append(drawdown_out, min(drawdown))
	    drawdown <- 0
            realized_profit <- append(realized_profit,
                                      (sum(coin_bag)*price) - sum(outflow)
            )
            # reset all values
            coin_bag <- c()
            buy_price <- c()
            outflow <- c()
            average_buy_price <- 0
            
            # buy ASAP
            coin_bag <- append(coin_bag, base_order_size/price)
            buy_prices <- append(buy_prices, price)
            outflow <- append(outflow, base_order_size)
            average_buy_price <- sum(outflow)/sum(coin_bag)
        }
        # sell all at the end of simulation
        if (i == len){
            TP_bool <- FALSE
    	    drawdown_out <- append(drawdown_out, min(drawdown))
            # update tibble for the last time
            headers <- 
                headers %>% 
                add_row(
                price = price,
                TP_bool = TP_bool,
                DCA_bool = DCA_bool,
                dca_orders_max = dca_orders_max,
                coin_bag = tail(coin_bag, n=1),
                buy_prices = tail(buy_prices, n=1),
                outflow = tail(outflow, n=1),
                average_buy_price = tail(average_buy_price, n=1),
                realized_profit = tail(realized_profit, n=1),
                price_deviation = price_deviation,
                dca_counter = dca_counter,
                price_deviation_target = price_deviation_target,
                step = i
            )
            # compute end of simulation results
#             profit_out <- headers %>% 
#                 filter(TP_bool == TRUE) %>% 
#                 pull(realized_profit) %>% 
#                 {sum(.)}
#             number_of_tp_actions <-
#                 headers %>% 
#                 filter(TP_bool == TRUE) %>% 
#                 pull(TP_bool) %>% 
#                 {length(.)}
# 	    uvc <- abs(mean(drawdown_out)) + (sd(drawdown_out)*1)
#             cost <- -(profit_out / (uvc/number_of_tp_actions))
        }
        # if (i %% 100 == 0){
        #     print(paste0("progress: ", round((i/len)*100,2), "%" ))
        #     }
    }
    # the lower the better
    return(headers)
}



list_of_runs <-
    runs %>% 
    group_by(iteration) %>% 
    group_split()

sample_runs <-
    list_of_runs %>% 
    {.[sample(c(1:length(list_of_runs)), size = n_simulations, replace = TRUE)]}



plan(multisession, workers = 8)
sum_cost_function <- function(
        theta
        ){
    out <-
	    # note: this runs are hardcoded, change this in the future
        sample_runs %>% 
        future_imap_dfr(
            ., function(X, idx){
            cost_sim <- cost_function(
                price_series = X,
                base_order_size = CAPITAL*theta[1],
                DCA_order_size = CAPITAL*theta[2],
                num_of_dca = 1,
                dca_multiplier = theta[3],
                price_deviation_target = theta[4],
                price_deviation_multiplier = theta[5],
                take_profit_target = theta[6],
                days = n_days,
                len = days
            )
            return(
                tibble(
                    cost = cost_sim,
                    idx = idx,
                    intercept = unique(X$intercept),
                    poly_1 = unique(X$poly_1),
                    poly_2 = unique(X$poly_2),
                    poly_3 = unique(X$poly_3),
                    mdl_sigma = unique(X$mdl_sigma)
                )
            )
            }
        )
    return(mean(out$cost))
}


# optimization ----


opts <- list(
    "algorithm" = "NLOPT_LN_BOBYQA",
    "xtol_rel" = 1.0e-04,
    "maxeval" = 5000000,
    "print_level" = 3
)

        # base_order_size,
        # DCA_order_size,
        # dca_multiplier,
        # num_of_dca,
        # price_deviation_target,
        # price_deviation_multiplier,
        # take_profit_target,
        # num_of_dca

print(paste0("Starting optimization for: ", coin))
non_linear_opt0 <-
    nloptr(
        x0 = c(
            0.01, # base order size
            0.01, # DCA order size
            1, # dca multiplier
            -0.01, # price_deviation_target
            1, # price deviation multiplier
            0.05 # take profit target
        ),
        eval_f = sum_cost_function,
        lb = c(0.01, 0.01, 1, -0.1, 1, 0.05),
        ub = c(0.25, 0.25, 3, -0.005, 3, 0.1),
        opts = opts
    )

optimization_out <-
    tibble(
        solution = non_linear_opt0$solution,
        parameters = c(
            "base_order_size",
            "DCA_order_size",
            "dca_multiplier",
            "price_deviation_target",
            "price_deviation_multiplier",
            "take_profit_target"
        )
    )

print(optimization_out)

write_csv(optimization_out, paste0(
    "optimization_results/optimization_",
    "coin_", coin,
    "days_", n_days,
    "capital_", CAPITAL,
    ".csv"))
#optimization_out <- read_csv(paste0("optimization_results/optimization_", coin, ".csv"))

# simulate with optimal parameters ----

sim_optimal <-
    sample_runs %>% 
    future_imap_dfr(
        ., function(X, idx){
        cost_sim <- cost_function2(
            price_series = X,
            base_order_size = CAPITAL*optimization_out$solution[1],
            DCA_order_size = CAPITAL*optimization_out$solution[2],
            num_of_dca = 1,
            dca_multiplier = optimization_out$solution[3],
            price_deviation_target = optimization_out$solution[4],
            price_deviation_multiplier = optimization_out$solution[5],
            take_profit_target = optimization_out$solution[6],
            days = n_days,
            len = n_days
        ) %>%
            mutate(
                idx = idx
                )
        }, .progress = TRUE)


profit <-
    sim_optimal %>% 
    ungroup() %>% 
    group_by(idx) %>% 
    filter(TP_bool == TRUE) %>% 
    summarise(
        percent_profit = (sum(realized_profit))/CAPITAL
    ) 
profit

write_csv(sim_optimal, paste0(
    "data/simulation_with_optimal_parameters_",
    "coin_", coin,
    "days_", n_days,
    "capital_", CAPITAL,
    ".csv"))

# print profit distribution
print(summary(profit$percent_profit))

# p1 <-
#     sim_optimal %>%
#     filter(price != 0, idx %in% c(1:20), average_buy_price != 0) %>%
#     ggplot(aes(
#         step, price, group = idx
#     )) +
#     geom_line() +
#     geom_line(aes(step, average_buy_price), color = "blue") +
#     facet_wrap(~idx)
# p1

    
