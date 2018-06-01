functions{
	#include "CIPlib.stan"
}

data{
	int T;
	real y[T];
	real sigma_x;
	real phi;
	real pri_mean;
}


parameters{
	real lam_y_bar;
	vector[T] x_bar;
}

transformed parameters{
	real lam_y;
	real prec_x;
	real sigma_y;
	real obsprec;
	
	vector[T] x;
	vector[T] y_scaled;
	vector[2*T] L;
	
	lam_y = lam_y_bar/sqrt(0.5*T);
	sigma_y = exp(-0.5*lam_y);
	prec_x = 1.0/pow(sigma_x,2);
	obsprec = 1.0/pow(sigma_y,2);
	for(i in 1:T){
		y_scaled[i] = obsprec*y[i];
	}
	L = CIP_TriDiagChol_const1n(T,prec_x + obsprec, prec_x*(1.0+pow(phi,2)) + obsprec,-phi*prec_x);
	x = CIP_TriDiagChol_LT_solve(L,x_bar) + CIP_TriDiagChol_LLT_solve(L,y_scaled);
	
	
	
}

model{
	target += normal_lpdf(lam_y | pri_mean, 3.0);
	target += normal_lpdf( x[1] | 0.0, sigma_x/sqrt(1.0-phi*phi));
	for( t in 2:T){
		target += normal_lpdf( x[t] | phi*x[t-1], sigma_x);
	}
	target += normal_lpdf(y | x, sigma_y);
	target += -L[2*T];
}


