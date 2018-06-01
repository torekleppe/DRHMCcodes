require(CIPlib)
require(rstan)
require(latex2exp)

y <- as.vector(read.table("USdata_updated.txt")$x)

T <- length(y)





cc <- stanc_builder(file="sw_drhmc_tauc.stan",allow_undefined=TRUE, isystem=CIP_header_path())
drhmc.mdlt <- stan_model(stanc_ret=cc,allow_undefined=TRUE,includes=CIP_include())


standta <- list(T=T,y=y,alpha=5.0,beta=0.5,zmean=0.0*y[1:(T-1)],xmean=0*y,lambdamean=0.0,lambdaprec=5+T-1.5)

opts <- optimizing(drhmc.mdlt,data=standta,as_vector=FALSE)
inits <- list(opts,opts,opts,opts,opts,opts,opts,opts,opts,opts)
drhmc <- sampling(drhmc.mdlt,data=standta,chains=10,seed=2,init=inits,control=list(stepsize_jitter=0.5))



summar <- (summary(drhmc,c("lambda","lambda_bar","z_bar[1]","x_bar[1]","tau_bar[1]","z[1]","x[1]","tau[1]"))$summary)

print(summar[,c(1,3,9)],digits=1)

print(mean(get_elapsed_time(drhmc)[,2]))

pdf("sw_plot_noprecomp.pdf",width=12,height=7)

time <- seq(from=1955+0.5*(0.25+0.5),to=2018+0.5*(0+0.25),length.out=T)

par(mfrow=c(1,3),ps=17)
#par(cex=1.1)
tauf <- summary(drhmc,"tau")$summary
plot(time,tauf[,"mean"],col=2,type="l",lwd=2,xlab="time (year)", ylab=TeX("$\\tau$"),
xlim=c(1955,2019),ylim=c(min(y),max(y)))
lines(time,tauf[,"2.5%"],col=4,lty=5)
lines(time,tauf[,"97.5%"],col=1,lty=4)
lines(time,y,type="p")
legend(x=1955,y=-1,c("mean","0.025-quantile","0.975-quantile","y")
,merge=TRUE,,col=c(2,4,1,1),lty=c(1,5,4,NA),lwd=c(2,1,1,NA),pch=c(NA,NA,NA,1))

xx <- summary(drhmc,"x")$summary
plot(time,xx[,"mean"],col=2,type="l",lwd=2,xlab="time (year)", ylab=TeX("$\\mathbf{x}$"),
xlim=c(1955,2019),ylim=c(-5,1.5))
lines(time,xx[,"2.5%"],col=4,lty=5)
lines(time,xx[,"97.5%"],col=1,lty=4)
legend(x=1955,y=1.5,c("mean","0.025-quantile","0.975-quantile")
,merge=TRUE,,col=c(2,4,1),lty=c(1,5,4),lwd=c(2,1,1))

zz <- summary(drhmc,"z")$summary
plot(time[1:T-1],zz[,"mean"],col=2,type="l",lwd=2,xlab="time (year)", ylab=TeX("$\\mathbf{z}$"),
xlim=c(1955,2019),ylim=c(-14,1))
lines(time[1:T-1],zz[,"2.5%"],col=4,lty=5)
lines(time[1:T-1],zz[,"97.5%"],col=1,lty=4)
legend(x=1955,y=-10,c("mean","0.025-quantile","0.975-quantile")
,merge=TRUE,,col=c(2,4,1),lty=c(1,5,4),lwd=c(2,1,1))
dev.off()

pdf("sw_trace_noprecomp.pdf",width=7,height=12)
par(mfrow=c(4,1))
chain <- extract(drhmc,c("lambda_bar"))$lambda_bar
TT <- length(chain)
plot(1:TT,chain,type="l",xlab="Iteration",ylab=TeX("$\\bar{ \\lambda}$"))
for(i in 0:10){
    lines(c(i*1000,i*1000),c(-5,5),col=2)
}


chain <- extract(drhmc,c("tau_bar[1]"))$"tau_bar[1]"
plot(1:TT,chain,type="l",xlab="Iteration",ylab=TeX("$\\bar{ \\tau}_1$"))
for(i in 0:10){
    lines(c(i*1000,i*1000),c(-5,5),col=2)
}


chain <- extract(drhmc,c("x_bar[1]"))$"x_bar[1]"
plot(1:TT,chain,type="l",xlab="Iteration",ylab=TeX("$\\bar{ \\mathbf{x}}_1$"))
for(i in 0:10){
    lines(c(i*1000,i*1000),c(-5,5),col=2)
}


chain <- extract(drhmc,c("z_bar[1]"))$"z_bar[1]"
plot(1:TT,chain,type="l",xlab="Iteration",ylab=TeX("$\\bar{ \\mathbf{z}}_1$"))
for(i in 0:10){
    lines(c(i*1000,i*1000),c(-5,5),col=2)
}


dev.off()



#stepsize_jitter=0.5,adapt_delta=0.8,metric="unit_e"
