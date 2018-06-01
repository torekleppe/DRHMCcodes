

data{
	int y[360];
	real fW[360];
	real mW[360];
	real WW[360];
	real fall[360];
	int f1[240];
	int m1[240];
	int f2[120];
	int m2[120];
}

transformed data{
	cov_matrix[2] fPri;
	fPri[1,1] = 1.0/1.244;
	fPri[1,2] = 0.0;
	fPri[2,1] = 0.0;
	fPri[2,2] = 1.0/1.244;
		

}

parameters{
	cov_matrix[2] WM;
	cov_matrix[2] WF;
	
	real<lower=0> kappaM;
	real<lower=0> kappaF;
	
	vector[2] b1f_raw[20];
	vector[2] b1m_raw[20];
	vector[20] b3f;
	vector[20] b3m;
	vector[5] beta;
}


transformed parameters{
	
	vector[40] b1f;
	vector[40] b1m;
	vector[360] eta;
	vector[2] zerovec;
	
	
	for( i in 1:20){
		b1f[i] = b1f_raw[i,1];
		b1f[i+20] = b1f_raw[i,2];
		b1m[i] = b1m_raw[i,1];
		b1m[i+20] = b1m_raw[i,2];
	}
	
	
	
	
	for(i in 1:240){
		eta[i] = beta[1] + beta[2]*fW[i] + beta[3]*mW[i] + beta[4]*WW[i] + beta[5]*fall[i] 
		+ b1f[f1[i]] + b1m[m1[i]]; 
	}
	
	for(i in 241:360){
		
		eta[i] = beta[1] + beta[2]*fW[i] + beta[3]*mW[i] + beta[4]*WW[i] + beta[5]*fall[i] 
		+ b3f[f2[i-240]] + b3m[m2[i-240]];
	}
	
	zerovec[1] = 0.0;
	zerovec[2] = 0.0;

}


model{
	target += wishart_lpdf(WM | 3.0,fPri);
	target += wishart_lpdf(WF | 3.0,fPri);
	
	
	target += multi_normal_prec_lpdf(b1f_raw | zerovec,WF);
	target += multi_normal_prec_lpdf(b1m_raw | zerovec,WM);
	
	
	
	target += gamma_lpdf(kappaM | 1.0, 0.622);
	target += gamma_lpdf(kappaF | 1.0, 0.622);
	
	target += normal_lpdf(b3f | 0.0, 1.0/sqrt(kappaF));
	target += normal_lpdf(b3m | 0.0, 1.0/sqrt(kappaM));

	target += bernoulli_logit_lpmf(y | eta);

}

generated quantities{
	real tau1f = (WF[1,1]*WF[2,2]-pow(WF[1,2],2))/WF[2,2];
	real tau1m = (WM[1,1]*WM[2,2]-pow(WM[1,2],2))/WM[2,2];
	real tau2f = (WF[1,1]*WF[2,2]-pow(WF[1,2],2))/WF[1,1];
	real tau2m = (WM[1,1]*WM[2,2]-pow(WM[1,2],2))/WM[1,1];
	real rhof = -WF[1,2]/sqrt(WF[1,1]*WF[2,2]);
	real rhom = -WM[1,2]/sqrt(WM[1,1]*WM[2,2]);
}





