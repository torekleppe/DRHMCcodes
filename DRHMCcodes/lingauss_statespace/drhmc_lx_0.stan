functions{
	#include "CIPlib.stan"
}

data{
	int T;
	real y[T];
	real sigma_y;
	real phi;
}

parameters{
	real lam_x_bar;
	vector[T] x_bar;
}

transformed parameters{
	real lam_x;
	real prec_x;
	real sigma_x;
	real obsprec;
	
	vector[T] x;
	vector[2*T] L;
	
	lam_x = lam_x_bar/sqrt(0.5*T);
	sigma_x = exp(-0.5*lam_x);
	prec_x = 1.0/pow(sigma_x,2);
	obsprec = 1.0/pow(sigma_y,2);
	
	L = CIP_TriDiagChol_const1n(T,prec_x + obsprec, prec_x*(1.0+pow(phi,2)) + obsprec,-phi*prec_x);
	x = CIP_TriDiagChol_LT_solve(L,x_bar);
	
	
}

model{

	target += normal_lpdf( x[1] | 0.0, sigma_x/sqrt(1.0-phi*phi));
	for( t in 2:T){
		target += normal_lpdf( x[t] | phi*x[t-1], sigma_x);
	}
	target += normal_lpdf(y | x, sigma_y);
	target += -L[2*T];
}


