
data{
	int T;
	real y[T];
}

parameters{
	
	real<lower=0> prec;
	real<lower=0,upper=1> phi_shifted;
	real mu;
	vector[T] zz;
}

transformed parameters{
	real sigma;
	real phi;
	
	vector[T] x;
	
	sigma = 1.0/sqrt(prec);
	phi = -1.0 + 2.0*phi_shifted;
	
	x[1] = mu + sigma/sqrt(1.0-phi*phi)*zz[1];
	for(t in 2:T){
		x[t] = mu + phi*(x[t-1]-mu) + sigma*zz[t];
	}
}

model{
	target += gamma_lpdf(prec | 5.0, 0.05);
	target += beta_lpdf(phi_shifted | 20.0, 1.5);
	target += normal_lpdf(mu | 0.0, 10.0);
	
	target += normal_lpdf(zz | 0.0, 1.0);
	target += normal_lpdf(y | 0.0, exp(0.5*x)); 

}