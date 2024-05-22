pacman::p_load(
    tidyverse,
    ggplot2
)

params <- read_csv("../optimization_results/optimization.csv")
params

order_size <- function(
        base_order_size,
        DCA_order_size,
        dca_multiplier,
        price_deviation_target,
        n_orders,
        capital
){
    c <- c()
    pd <- c(1, cumsum(abs(rep(price_deviation_target, n_orders)))[1:n_orders-1])
    quantity <- c()
    average_buy_price <- c()
    for (i in 1:n_orders){
        if (i < 3){
            dm <- 1
        }
        else{
            dm <- dca_multiplier
        }
        calc <- (capital*base_order_size) * as.numeric(i==1) +
            (( (capital*DCA_order_size) * (dm^(i-2))) * as.numeric(i>1))
        c <- append(c, calc)
        quantity <- append(quantity, calc/pd[i])
        average_buy_price <- append(average_buy_price, sum(quantity*pd[1:i])/sum(quantity))
    }
    return(pd)
}

order_size(
    0.010173,
    0.01007,
    1.1,
    -0.0133,
    100,
    23259.66661347
)
%>% plot(., type = "l")
