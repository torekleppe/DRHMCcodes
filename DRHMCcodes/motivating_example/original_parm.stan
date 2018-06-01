data{
	real y;
}
parameters{
	real q1;
	real q2;
}
model{
	// "priors"
	target += normal_lpdf(q1 | 0.0 , 1.0);
	target += normal_lpdf(q2 | 0.0 , 1.0);
	// "likelihood"
	target += normal_lpdf(y | q2, exp(-1.5*q1)); // notice: standard deviation
}

