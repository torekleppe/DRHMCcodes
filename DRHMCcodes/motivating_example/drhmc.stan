data{
	real y;
}

/* 
	Change sampled parameters to their standardised/bar counterparts in the
	parameters block
*/
parameters{
	real q1_bar;
	real q2_bar;
}

/*
	Transformation from modified to original parameters is  
	implemented in the (new) "transformed parameters"-block
*/
transformed parameters{
	real q1;
	real q2;
	real h2;
	real L2;
	// first block variable transformation: G_(1)=1+9/2
	q1 = q1_bar/sqrt(1.0+4.5);
	
	// now q1 is available, can compute h_(2) and L_(2)
	h2 = y/(1.0+exp(-3.0*q1)); // = E(q2|q1,y)
	L2 = sqrt(1.0+exp(3.0*q1));
	
	// second block variable transformation:
	q2 = h2 + q2_bar/L2;
}

/* 	
	Model specification is identical to the original code,
	except for contribution from Jacobian 
*/
model{ 
	// "priors"
	target += normal_lpdf(q1 | 0.0 , 1.0);
	target += normal_lpdf(q2 | 0.0 , 1.0);
	// "likelihood"
	target += normal_lpdf(y | q2, exp(-1.5*q1)); // notice: standard deviation
	// contribution from Jacobian
	target += -log(L2);
}

