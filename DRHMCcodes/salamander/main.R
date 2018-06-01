require(rstan)

## load the salamander data
load("salam.RData")

## organize data into a form suitable for logistic regression
dat0=data.frame("y"=c(salam$y),
"fW"=as.integer(salam$x[,"W/R"]==1 | salam$x[,"W/W"]==1),
"mW"=as.integer(salam$x[,"R/W"]==1 | salam$x[,"W/W"]==1),
"WW"=as.integer(salam$x[,"W/W"]==1 ) )

## add salamander id (female, male)
id = t( apply(salam$z, 1, function(x) {
    tmp = which (x==1)
    tmp[2] = tmp[2] - 20
    tmp
}) )

## index for the experiment group
dat0$group=as.factor(c(rep(c(rep(1,5),rep(2,10),rep(1,5)),6),
rep(c(rep(3,5),rep(4,10),rep(3,5)),6),
rep(c(rep(5,5),rep(6,10),rep(5,5)),6)))
## index for the experiment
dat0$experiment=as.factor(rep(1:3, each=120))

## set the indices for the first two experiments modeled as iid2d,
## (The indices for the third experiment are set to NA)
fid_iid2_e1e2 <- c(id[,1],id[,1] + 20, rep(NA, 120))
mid_iid2_e1e2 <- c(id[,2],id[,2] + 20, rep(NA, 120))

## set the indices for third experiment
## (The indices for the first two experiments are set to NA)
fid_id_e3 <- c(rep(NA,240),id[,1])
mid_id_e3 <- c(rep(NA,240),id[,2])

## Indicator for fall
fall <- c(rep(0, 120), rep(1,240))

## generate the dataset
data <- data.frame(dat0, fid_iid2_e1e2, mid_iid2_e1e2, fid_id_e3, mid_id_e3, fall)

standta.0 <- list(y=data$y,fW=data$fW,mW=data$mW,WW=data$WW,fall=data$fall,f1=(data$fid_iid2_e1e2)[1:240],
m1=(data$mid_iid2_e1e2)[1:240],f2=(data$fid_id_e3)[241:360],m2=(data$mid_id_e3)[241:360]);

#default<-stan(file="salamander_default.stan",data=standta.0,chains=10,control=list(adapt_delta=0.95),seed=1)


standta <- list(y=data$y,fW=data$fW,mW=data$mW,WW=data$WW,fall=data$fall,f1=(data$fid_iid2_e1e2)[1:240],
m1=(data$mid_iid2_e1e2)[1:240],f2=(data$fid_id_e3)[241:360],m2=(data$mid_id_e3)[241:360],
lamFmean=c(0.0,0.0),lamMmean=c(0.0,0.0),V1Fmean=0.0,V1Mmean=0.0,logkappaFmean=0.0,logkappaMmean=0.0);

drhmc.latents<-stan(file="salamander_drhmc_latents.stan",data=standta,chains=10,control=list(adapt_delta=0.95),seed=1)


drhmc.wish<-stan(file="salamander_drhmc_wish.stan",data=standta,chains=10,control=list(adapt_delta=0.95),seed=1)


tab.def <-(summary(default,c("tau1f","tau2f","rhof","tau1m","tau2m","rhom","kappaF","kappaM","beta"))$summary)
tab.lat <-(summary(drhmc.latents,c("tau1f","tau2f","rhof","tau1m","tau2m","rhom","kappaF","kappaM","beta"))$summary)
tab.wish <-(summary(drhmc.wish,c("tau1f","tau2f","rhof","tau1m","tau2m","rhom","kappaF","kappaM","beta"))$summary)


tab <- cbind(tab.def[,c(1,3,9)],tab.lat[,c(1,3,9)],tab.wish[,c(1,3,9)])

print(tab,digits=1)









#lamFmean <- (summary(drhmc,c("lamF"))$summary)[,1];
#lamMmean <- (summary(drhmc,c("lamM"))$summary)[,1]
#V1Fmean <- (summary(drhmc,c("V1F"))$summary)[,1]
#V1Mmean <- (summary(drhmc,c("V1M"))$summary)[,1]
#logkappaFmean <- (summary(drhmc,c("logkappaF"))$summary)[,1]
#logkappaMmean <- (summary(drhmc,c("logkappaM"))$summary)[,1]

#standta <- list(y=data$y,fW=data$fW,mW=data$mW,WW=data$WW,fall=data$fall,f1=(data$fid_iid2_e1e2)[1:240],
#m1=(data$mid_iid2_e1e2)[1:240],f2=(data$fid_id_e3)[241:360],m2=(data$mid_id_e3)[241:360],
#lamFmean=lamFmean,lamMmean=lamMmean,V1Fmean=V1Fmean,V1Mmean=V1Mmean,logkappaFmean=logkappaFmean,
#logkappaMmean=logkappaMmean);

#drhmc.mean<-stan(file="salamander_drhmc_wish.stan",data=standta,chains=10,control=list(adapt_delta=0.95),seed=1)
