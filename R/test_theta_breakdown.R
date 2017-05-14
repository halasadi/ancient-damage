

####################  test theta breakdown  ################################


model <- get(load("../utilities/modern_Jones/clus_2/model.rda"))
theta_pool <- model$theta
ll <- theta_breakdown(theta_pool)


#################  test theta breakdown Lazaridis + Skoglund + moderns  ##################


model <- get(load("../utilities/modern_Skoglund_Lazaridis/clus_3/model.rda"))
theta_pool <- model$theta
ll <- theta_breakdown(theta_pool)

a <- ll[[2]]
save(a, file = "../utilities/logo_design_example_1.rda")


#################  test theta breakdown modern + Pinhasi + Lazaridis ##################


model <- get(load("../utilities/modern_Pinhasi_Lazaridis/clus_3/model.rda"))
theta_pool <- model$theta
ll <- theta_breakdown(theta_pool)

b <- ll[[1]]
save(b, file = "../utilities/logo_design_example_2.rda")

