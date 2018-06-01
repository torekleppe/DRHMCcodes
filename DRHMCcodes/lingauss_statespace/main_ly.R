require(CIPlib)
require(rstan)
require(latex2exp)


# high observation noise case
lam_x <- -log(0.15^2)
lam_y <- -log(0.005^2)
omega <- 2.2

source("sim_data.R")
source("kalman_filters.R")

set1 <- read.table("high_obs_noise.txt")
T <- length(set1$y)


#
# lam_x experiment
#

kern <- post_kern_ly(from=0,to=25,ng=1000,lam_x=lam_x,omega=omega,set1$y, pri_mean=lam_y)

standta <- list(T = T, y = set1$y, sigma_x = sigma_x, phi = phi, pri_mean = lam_y )




cout <- stanc_builder(file="drhmc_ly_0.stan",allow_undefined=TRUE,isystem=CIP_header_path())
smodel <- stan_model(stanc_ret=cout,allow_undefined=TRUE,include=CIP_include())
drhmc.0 <- sampling(smodel,chains=10,seed=1,data=standta)
cout <- stanc_builder(file="drhmc_ly_y.stan",allow_undefined=TRUE,isystem=CIP_header_path())
smodel <- stan_model(stanc_ret=cout,allow_undefined=TRUE,include=CIP_include())
drhmc.y <- sampling(smodel,chains=10,seed=1,data=standta)
cout <- stanc_builder(file="drhmc_ly_xcy.stan",allow_undefined=TRUE,isystem=CIP_header_path())
smodel <- stan_model(stanc_ret=cout,allow_undefined=TRUE,include=CIP_include())
drhmc.xcy <- sampling(smodel,chains=10,seed=1,data=standta)
hmc <- stan(file="hmc_ly.stan",data=standta,chains=10,seed=1)


pdf("ly_plot_low.pdf",width=14,height=7)
par(mfrow=c(2,4))
ts.plot(extract(hmc,c("lam_y"))$lam_y,ylab=TeX("$\\tau$"),xlab="Iteration",main=TeX("$\\mathbf{x}$-prior standardisation"),xlim=c(0,10000),ylim=c(5.0,25.0))
for(i in 0:10){
    lines(c(i*1000,i*1000),c(0,100),col=2)
}
ts.plot(extract(drhmc.0,c("lam_y"))$lam_y,ylab=TeX("$\\tau$"),xlab="Iteration",
main=TeX("DRHMC, $\\mathbf{h}_{(2)} = 0$"),xlim=c(0,10000),ylim=c(5.0,25.0))
for(i in 0:10){
    lines(c(i*1000,i*1000),c(0,100),col=2)
}
ts.plot(extract(drhmc.y,c("lam_y"))$lam_y,ylab=TeX("$\\tau$"),xlab="Iteration",main=TeX("DRHMC, $\\mathbf{h}_{(2)} = \\mathbf{y}$"),xlim=c(0,10000),ylim=c(5.0,25.0))
for(i in 0:10){
    lines(c(i*1000,i*1000),c(0,100),col=2)
}

ts.plot(extract(drhmc.xcy,c("lam_y"))$lam_y,ylab=TeX("$\\tau$"),xlab="Iteration",
main=TeX("DRHMC, $\\mathbf{h}_{(2)}=E(\\mathbf{x}|\\mathbf{y},\\lambda,\\tau)$"),xlim=c(0,10000),ylim=c(5.0,25.0))
for(i in 0:10){
    lines(c(i*1000,i*1000),c(0,100),col=2)
}


hist(extract(hmc,c("lam_y"))$lam_y,100,probability=T,main="",xlab=TeX("$\\tau$"), xlim=c(5.0,25.0))
lines(x=kern[,1],y=kern[,2]/(kern[2,1]-kern[1,1]),col=2)
hist(extract(drhmc.0,c("lam_y"))$lam_y,100,probability=T,main="",xlab=TeX("$\\tau$"), xlim=c(5.0,25.0))
lines(x=kern[,1],y=kern[,2]/(kern[2,1]-kern[1,1]),col=2)
hist(extract(drhmc.y,c("lam_y"))$lam_y,100,probability=T,main="",xlab=TeX("$\\tau$"), xlim=c(5.0,25.0))
lines(x=kern[,1],y=kern[,2]/(kern[2,1]-kern[1,1]),col=2)
hist(extract(drhmc.xcy,c("lam_y"))$lam_y,100,probability=T,main="",xlab=TeX("$\\tau$"), xlim=c(5.0,25.0))
lines(x=kern[,1],y=kern[,2]/(kern[2,1]-kern[1,1]),col=2)
dev.off()

true.mean <- sum(kern[,1]*kern[,2])
true.m2 <- sum(kern[,1]^2*kern[,2])
true.std <- sqrt(true.m2-true.mean^2)

tab <- matrix(0.0,nrow=5,ncol=6)
tab[1,1] <- true.mean
tab[1,2] <- true.std
tab[2,1:3] <- (summary(hmc,c("lam_y"))$summary)[c(1,3,9)]
tab[2,4] <- mean(get_elapsed_time(hmc)[,2])
tab[3,1:3] <- (summary(drhmc.0,c("lam_y"))$summary)[c(1,3,9)]
tab[3,4] <- mean(get_elapsed_time(drhmc.0)[,2])
tab[4,1:3] <- (summary(drhmc.y,c("lam_y"))$summary)[c(1,3,9)]
tab[4,4] <- mean(get_elapsed_time(drhmc.y)[,2])
tab[5,1:3] <- (summary(drhmc.xcy,c("lam_y"))$summary)[c(1,3,9)]
tab[5,4] <- mean(get_elapsed_time(drhmc.xcy)[,2])




