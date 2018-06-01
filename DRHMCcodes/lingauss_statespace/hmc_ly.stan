data{
	int T;
	real y[T];
	real sigma_x;
	real phi;
	real pri_mean;
}

parameters{
	real lam_y;
	vector[T] zz;
}

transformed parameters{
	
	real sigma_y;
	vector[T] x;
	
	sigma_y = exp(-0.5*lam_y);
	x[1] = sigma_x/sqrt(1.0-phi*phi)*zz[1];
	for( t in 2:T){
		x[t] = phi*x[t-1] + sigma_x*zz[t];
	}
	
	
}

model{
	target += normal_lpdf(lam_y | pri_mean, 3.0);
	target += normal_lpdf(zz | 0.0, 1.0);
	target += normal_lpdf(y | x, sigma_y);
}


