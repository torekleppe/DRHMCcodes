data{
	int T;
	real y[T];
	real sigma_y;
	real phi;
}

parameters{
	real lam_x;
	vector[T] zz;
}

transformed parameters{
	
	real sigma_x;
	vector[T] x;
	
	sigma_x = exp(-0.5*lam_x);
	x[1] = sigma_x/sqrt(1.0-phi*phi)*zz[1];
	for( t in 2:T){
		x[t] = phi*x[t-1] + sigma_x*zz[t];
	}
	
	
}

model{

	target += normal_lpdf( zz | 0.0, 1.0);
	target += normal_lpdf(y | x, sigma_y);
}


