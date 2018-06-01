data{
	int T;
	real y[T];
	real alpha;
	real beta;
}



parameters{
	real lambda;
	vector[T-1] z;
	vector[T] x;
	vector[T] tau;
}

transformed parameters{
	real sigma;
	sigma = exp(-0.5*lambda);

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
	
	
}


