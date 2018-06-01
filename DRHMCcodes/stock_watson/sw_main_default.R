
require(rstan)


y <- as.vector(read.table("USdata_updated.txt")$x)

T <- length(y)


standta <- list(T=T,y=y,alpha=5.0,beta=0.5)
default <- stan(file="sw_default.stan",data=standta,chains=10,seed=1,control=list(adapt_delta=0.99))

print(summary(default,c("lambda","z[1]","x[1]","tau[1]"))$summary)

print(mean(get_elapsed_time(default)[,2]))

