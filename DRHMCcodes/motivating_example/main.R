require(rstan)

data <- list(y=0.5)
original <- stan(file="original_parm.stan",data=data,control=list(adapt_delta=0.99))
drhmc <- stan(file="drhmc.stan",data=data,control=list(adapt_delta=0.9))

