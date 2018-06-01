set.seed(1)
T <- 100
x <- vector(length=T)
y <- x



sigma_x <- exp(-0.5*lam_x)
sigma_y <- exp(-0.5*lam_y)
phi <- tanh(CIP_AR1_psi(omega,T)$psi)

x[1] <- sigma_x/sqrt(1-phi*phi)*rnorm(1)
y[1] <- x[1]+sigma_y*rnorm(1)

for( t in 2:T){
    x[t] <- phi*x[t-1] + sigma_x*rnorm(1)
    y[t] <- x[t] + sigma_y*rnorm(1)
}

dta <- cbind(x,y)
write.table(x=dta,file="high_obs_noise.txt")

