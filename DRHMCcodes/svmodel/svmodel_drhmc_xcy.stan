functions{
	#include "CIPlib.stan"
}

data{
	int T;
	real y[T];
}

transformed data{
	real omega_defaultPrior_prec;
	vector[T] halfLogYsq;
	omega_defaultPrior_prec = CIP_AR1_omega_defaultPrior_prec(20.0,1.5,T);
	
	for(t in 1:T){
		halfLogYsq[t] = 0.5*log(pow(y[t],2));
	}
}

parameters{
	real lambda_bar;
	real omega_bar;
	real mu_bar;
	vector[T] x_bar;
}

transformed parameters{
	real lambda;
	real sigma;
	real precx;
	real omega;
	real psi;
	real phi;
	real L3;
	real mu;
	real tmp0;
	
	vector[2*T] L4;
	vector[T] x;
	vector[T] tmp;
	
	// block 1
	lambda = lambda_bar/sqrt(5.0 + 0.5*T);
	sigma = exp(-0.5*lambda);
	precx = exp(lambda);
	
	// block 2
	omega = omega_bar/sqrt(omega_defaultPrior_prec + 0.5*T);
	psi = CIP_AR1_psi(omega,T);
	phi = tanh(psi);
	
	// block 3
	L3 = sqrt(0.01 + exp(lambda)*(2.0*(T-1)*(1.0-phi)-(T-2)/pow(cosh(psi),2)));
	mu = mu_bar/L3;
	
	// block 4
	L4 =  CIP_TriDiagChol_const1n(T,precx + 0.5, precx*(1.0+pow(phi,2)) + 0.5,-phi*precx);
	tmp[1] = (1.0-phi)*mu*precx;
	tmp[T] = tmp[1];
	tmp0 = pow(1.0-phi,2)*mu*precx;
	for( t in 2:(T-1)){
		tmp[t] = tmp0;
	}
	tmp += halfLogYsq;
	x = CIP_TriDiagChol_LT_solve(L4,x_bar) + CIP_TriDiagChol_LLT_solve(L4,tmp);

}

model{
	target += CIP_ExpGammaR_logpdf(lambda,5.0,0.05);
	target += CIP_AR1_omega_defaultPrior_logpdf(omega,20.0,1.5,T);
	target += normal_lpdf(mu | 0.0, 10.0);
	
	// latent process
	// t=1
	target += normal_lpdf(x[1] | mu , sigma/sqrt(1.0-pow(phi,2)));
	
	// remaining times
	for( t in 2:T){
		target += normal_lpdf(x[t] | mu + phi*(x[t-1]-mu),sigma);
	}
	
	// observations
	
	target += normal_lpdf(y | 0.0,exp(0.5*x));
	
	// Jaobian correction

	target += -log(L3);
	target += -L4[2*T];
}


