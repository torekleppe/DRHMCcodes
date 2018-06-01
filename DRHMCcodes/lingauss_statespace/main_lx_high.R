require(CIPlib)
require(rstan)
require(latex2exp)


# high observation noise case
lam_x <- -log(0.15^2)
lam_y <- -log(0.15^2)
omega <- 2.2

source("sim_data.R")
source("kalman_filters.R")

set1 <- read.table("high_obs_noise.txt")
T <- length(set1$y)


#
# lam_x experiment
#

kern <- post_kern_lx(from=0,to=10,ng=1000,lam_y=lam_y,omega=omega,set1$y)

standta <- list(T = T, y = set1$y, sigma_y = sigma_y, phi = phi)




cout <- stanc_builder(file="drhmc_lx_0.stan",allow_undefined=TRUE,isystem=CIP_header_path())
smodel <- stan_model(stanc_ret=cout,allow_undefined=TRUE,include=CIP_include())
drhmc.0 <- sampling(smodel,chains=10,seed=1,data=standta)
cout <- stanc_builder(file="drhmc_lx_y.stan",allow_undefined=TRUE,isystem=CIP_header_path())
smodel <- stan_model(stanc_ret=cout,allow_undefined=TRUE,include=CIP_include())
drhmc.y <- sampling(smodel,chains=10,seed=1,data=standta)
cout <- stanc_builder(file="drhmc_lx_xcy.stan",allow_undefined=TRUE,isystem=CIP_header_path())
smodel <- stan_model(stanc_ret=cout,allow_undefined=TRUE,include=CIP_include())
drhmc.xcy <- sampling(smodel,chains=10,seed=1,data=standta)
hmc <- stan(file="hmc_lx.stan",data=standta,chains=10,seed=1)


pdf("lx_plot_high.pdf",width=14,height=7)
par(mfrow=c(2,4))
ts.plot(extract(hmc,c("lam_x"))$lam_x,ylab=TeX("$\\lambda_x$"),xlab="Iteration",main=TeX("prior rescaling"),xlim=c(0,10000),ylim=c(3.0,5.2))
for(i in 0:10){
    lines(c(i*1000,i*1000),c(0,100),col=2)
}
ts.plot(extract(drhmc.0,c("lam_x"))$lam_x,ylab=TeX("$\\lambda_x$"),xlab="Iteration",
main=TeX("DRHMC, $\\mu_{(2)} = 0$"),xlim=c(0,10000),ylim=c(3.0,5.2))
for(i in 0:10){
    lines(c(i*1000,i*1000),c(0,100),col=2)
}
ts.plot(extract(drhmc.y,c("lam_x"))$lam_x,ylab=TeX("$\\lambda_x$"),xlab="Iteration",main=TeX("DRHMC, $\\mu_{(2)} = \\mathbf{y}$"),xlim=c(0,10000),ylim=c(3.0,5.2))
for(i in 0:10){
    lines(c(i*1000,i*1000),c(0,100),col=2)
}

ts.plot(extract(drhmc.xcy,c("lam_x"))$lam_x,ylab=TeX("$\\lambda_x$"),xlab="Iteration",
main=TeX("DRHMC, $\\mu_{(2)}=E(\\mathbf{x}|\\mathbf{y},\\lambda_x)$"),xlim=c(0,10000),ylim=c(3.0,5.2))
for(i in 0:10){
    lines(c(i*1000,i*1000),c(0,100),col=2)
}


hist(extract(hmc,c("lam_x"))$lam_x,100,probability=T,main="",xlab=TeX("$\\lambda_x$"), xlim=c(3.0,5.2))
lines(x=kern[,1],y=kern[,2]/(kern[2,1]-kern[1,1]),col=2)
hist(extract(drhmc.0,c("lam_x"))$lam_x,100,probability=T,main="",xlab=TeX("$\\lambda_x$"), xlim=c(3.0,5.2))
lines(x=kern[,1],y=kern[,2]/(kern[2,1]-kern[1,1]),col=2)
hist(extract(drhmc.y,c("lam_x"))$lam_x,100,probability=T,main="",xlab=TeX("$\\lambda_x$"), xlim=c(3.0,5.2))
lines(x=kern[,1],y=kern[,2]/(kern[2,1]-kern[1,1]),col=2)
hist(extract(drhmc.xcy,c("lam_x"))$lam_x,100,probability=T,main="",xlab=TeX("$\\lambda_x$"), xlim=c(3.0,5.2))
lines(x=kern[,1],y=kern[,2]/(kern[2,1]-kern[1,1]),col=2)
dev.off()

true.mean <- sum(kern[,1]*kern[,2])
true.m2 <- sum(kern[,1]^2*kern[,2])
true.std <- sqrt(true.m2-true.mean^2)

tab <- matrix(0.0,nrow=5,ncol=6)
tab[1,1] <- true.mean
tab[1,2] <- true.std
tab[2,1:3] <- (summary(hmc,c("lam_x"))$summary)[c(1,3,9)]
tab[2,4] <- mean(get_elapsed_time(hmc)[,2])
tab[3,1:3] <- (summary(drhmc.0,c("lam_x"))$summary)[c(1,3,9)]
tab[3,4] <- mean(get_elapsed_time(drhmc.0)[,2])
tab[4,1:3] <- (summary(drhmc.y,c("lam_x"))$summary)[c(1,3,9)]
tab[4,4] <- mean(get_elapsed_time(drhmc.y)[,2])
tab[5,1:3] <- (summary(drhmc.xcy,c("lam_x"))$summary)[c(1,3,9)]
tab[5,4] <- mean(get_elapsed_time(drhmc.xcy)[,2])




