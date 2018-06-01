require(rstan)
require(CIPlib)

y <- read.table("returns_sp500_99_09.txt")
y <- 100*y$V1
T <- length(y)

standta <- list(T=T,y=y)



hmc <- stan(file="svmodel_default.stan",data=standta,chains=10,control=list(adapt_kappa=0.95,adapt_delta=0.9),seed=1)
table <- matrix(0.0,nrow=16,ncol=7)
table[1:3,1:5] <- t((summary(hmc,c("sigma","phi","mu","x[1]","x[2515]"))$summary)[,c(1,3,9)])
table[1,6] <- mean(get_elapsed_time(hmc)[,2])
table[4,1:5] <- table[3,1:5]/sum(get_elapsed_time(hmc)[,2])
print(warnings())
rm(hmc)


cc <- stanc_builder(file="svmodel_drhmc_mu.stan",allow_undefined=TRUE,isystem=CIP_header_path())
mdl.mu <- stan_model(stanc_ret=cc,allow_undefined=TRUE, include=CIP_include())
drhmc.mu <- sampling(mdl.mu,data=standta,chains=10,control=list(adapt_kappa=0.95,adapt_delta=0.9),seed=1)

table[5:7,1:5] <- t((summary(drhmc.mu,c("sigma","phi","mu","x[1]","x[2515]"))$summary)[,c(1,3,9)])
table[5,6] <- mean(get_elapsed_time(drhmc.mu)[,2])
table[8,1:5] <- table[7,1:5]/sum(get_elapsed_time(drhmc.mu)[,2])
rm(drhmc.mu)
print(warnings())

cc <- stanc_builder(file="svmodel_drhmc_xcy.stan",allow_undefined=TRUE,isystem=CIP_header_path())
mdl.xcy <- stan_model(stanc_ret=cc,allow_undefined=TRUE, include=CIP_include())
drhmc.xcy <- sampling(mdl.xcy,data=standta,chains=10,control=list(adapt_kappa=0.95,adapt_delta=0.9),seed=1)

table[9:11,1:5] <- t((summary(drhmc.xcy,c("sigma","phi","mu","x[1]","x[2515]"))$summary)[,c(1,3,9)])
table[9,6] <- mean(get_elapsed_time(drhmc.xcy)[,2])
table[12,1:5] <- table[11,1:5]/sum(get_elapsed_time(drhmc.xcy)[,2])

print(warnings())

cc <- stanc_builder(file="svmodel_drhmc_0.stan",allow_undefined=TRUE,isystem=CIP_header_path())
mdl.0 <- stan_model(stanc_ret=cc,allow_undefined=TRUE, include=CIP_include())
drhmc.0 <- sampling(mdl.0,data=standta,chains=10,control=list(adapt_kappa=0.95,adapt_delta=0.9),seed=1)

table[13:15,1:5] <- t((summary(drhmc.0,c("sigma","phi","mu","x[1]","x[2515]"))$summary)[,c(1,3,9)])
table[13,6] <- mean(get_elapsed_time(drhmc.0)[,2])
table[16,1:5] <- table[15,1:5]/sum(get_elapsed_time(drhmc.0)[,2])
print(warnings())


write.table(table,file="table.txt")

#control=list(adapt_kappa=0.95,adapt_delta=0.9)
