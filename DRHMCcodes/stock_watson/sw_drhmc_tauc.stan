functions{
	#include "CIPlib.stan"
}

data{
	int T;
	real y[T];
	real alpha;
	real beta;
	real lambdamean;
	real lambdaprec;
	real zmean[T-1];
	real xmean[T];
}



parameters{
	real lambda_bar;
	vector[T-1] z_bar;
	vector[T] x_bar;
	vector[T] tau_bar;
}

transformed parameters{
	real lambda;
	real sigma;
	real lprec;
	vector[T-1] z;
	vector[T] x;
	vector[T] tau;
	
	vector[2*(T-1)] Lz;
	vector[2*T] Lx;
	vector[T-1] expmz;
	vector[T] G4diag;
	vector[2*T] Ltau;
	vector[T] Pyy;
	
	// block 1
	lambda =  lambdamean + lambda_bar/sqrt(lambdaprec);
	sigma = exp(-0.5*lambda);
	lprec = exp(lambda);
	
	// block 2
	Lz = CIP_TriDiagChol_const1n(T-1,lprec+0.5,2.0*lprec+0.5,-lprec);
	z = CIP_TriDiagChol_LT_solve(Lz,z_bar);
	for(t in 1:(T-1)){
		z[t] += zmean[t];
	}
	
	// block 3
	Lx = CIP_TriDiagChol_const1n(T,lprec+0.5,2.0*lprec+0.5,-lprec);
	x = CIP_TriDiagChol_LT_solve(Lx,x_bar);
	for(t in 1:T){
		x[t] += xmean[t];
	}
	
	// block 4
	expmz = exp(-z);
	G4diag = exp(-x);
	for(t in 1:T){
		Pyy[t] = G4diag[t]*y[t];
	}
	
	G4diag[1] += expmz[1];
	for(t in 2:(T-1)){
		G4diag[t] += expmz[t-1]+expmz[t];
	}
	G4diag[T] += expmz[T-1];
	
	
	
	Ltau = CIP_TriDiagChol(G4diag,-expmz);
	tau = CIP_TriDiagChol_LT_solve(Ltau,tau_bar) + CIP_TriDiagChol_LLT_solve(Ltau,Pyy);
	
	
	
	
			 
}


model{
	target += alpha*lambda - beta*exp(lambda);
	
	for( t in 2:(T-1)) {
		target += normal_lpdf(z[t] | z[t-1], sigma);
	}
	
	for( t in 2:T){
		target += normal_lpdf(x[t] | x[t-1], sigma);
	}
	
	for( t in 2:T){
		target += normal_lpdf(tau[t] | tau[t-1],exp(0.5*z[t-1]));
	}
	
	for( t in 1:T){
		target += normal_lpdf(y[t] | tau[t], exp(0.5*x[t]));
	}
	
	target += -Lz[2*(T-1)] -Lx[2*T] -Ltau[2*T];
	
}




