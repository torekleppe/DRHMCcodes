require(CIPlib)
require(stats)
posterior_kernel <- function(lam_x,lam_y,omega,y){
    
    T <- length(y);
    phi <- tanh(CIP_AR1_psi(omega,T)$psi)
    sigmav <- exp(-0.5*lam_x);
    sigmay <- exp(-0.5*lam_y);
    margstd <- sigmav/sqrt(1.0-phi^2);
    
    # kalman filter
    aa <- 0.0;
    PP <- margstd^2;
    ll <- 0.0
    for (t in 1:T){
        et <- y[t] - aa;
        Dt <- PP + sigmay^2
        ll <- ll - 0.5*log(Dt) - 0.5*et^2/Dt;
        Kt <- (phi*PP)/Dt;
        aa <- phi*aa + Kt*et;
        Lt <- phi - Kt;
        Jt <- sigmav; # - Kt*sigmay;
        PP <- phi*PP*Lt + sigmav*Jt;
    }
    return(ll)
}



post_kern_lx <- function(from,to,ng,lam_y,omega,y){
    grid <- as.vector(seq(from=from,to=to,length.out=ng))
    lkern <- 0*grid
    mu.1 <- 0*grid
    mu.T <- 0*grid
    sig.1 <- 0*grid
    sig.T <- 0*grid
    
    T <- length(y);
    phi <- tanh(CIP_AR1_psi(omega,T)$psi)
    sigmay <- exp(-0.5*lam_y);
    
    for(i in 1:ng){
        
        lkern[i] <- posterior_kernel(grid[i],lam_y,omega,y)
        #sigmax <- exp(-0.5*grid[i]);
        #Sig <- toeplitz(sigmax^2/(1.0-phi^2)*(phi^(0:(T-1))))
        
        #mu.post <- Sig%*%solve(Sig+sigmay^2*diag(1.0,T,T),y)
        #Sig.post <- Sig - Sig%*%solve(Sig+sigmay^2*diag(1.0,T,T),Sig)
        #mu.1[i] <- mu.post[1]
        #mu.T[i] <- mu.post[T]
        #sig.1[i] <- sqrt(Sig.post[1,1])
        #sig.T[i] <- sqrt(Sig.post[T,T])
        
    }
    
    wts <- exp(lkern-max(lkern))
    wts <- wts/sum(wts)
    
    ret <- cbind(grid,wts,lkern)#,mu.1,sig.1,mu.T,sig.T)
    
}

post_kern_ly <- function(from,to,ng,lam_x,omega,y,pri_mean){
    grid <- as.vector(seq(from=from,to=to,length.out=ng))
    lkern <- 0*grid
    
    
    T <- length(y);
    phi <- tanh(CIP_AR1_psi(omega,T)$psi)
    
    for(i in 1:ng){
        
        lkern[i] <- posterior_kernel(lam_x,grid[i],omega,y)  - 0.5*(grid[i]-pri_mean)^2/(9.0)
        
    }
    
    wts <- exp(lkern-max(lkern))
    wts <- wts/sum(wts)
    
    ret <- cbind(grid,wts,lkern)
    
}







