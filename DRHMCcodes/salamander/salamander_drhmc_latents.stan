

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
	
	vector[40] zf1;
	vector[40] zm1;
	vector[20] zf3;
	vector[20] zm3;
	vector[5] beta;
}


transformed parameters{
	vector[40] b1f;
	vector[40] b1m;
	vector[20] b3f;
	vector[20] b3m;
	
	vector[360] eta;
	
	real detf;
	real wt1f;
	real wt12f;
	real wt2f;
	real detm;
	real wt1m;
	real wt12m;
	real wt2m;
	
	real stdkm;
	real stdkf;
	
	
	
	detf = WF[1,1]*WF[2,2] - pow(WF[1,2],2);
	wt1f = sqrt(WF[2,2]/detf);
	wt12f = - WF[1,2]/(detf*wt1f);
	wt2f = sqrt(1.0/WF[2,2]);
	
	detm = WM[1,1]*WM[2,2] - pow(WM[1,2],2);
	wt1m = sqrt(WM[2,2]/detm);
	wt12m = - WM[1,2]/(detm*wt1m);
	wt2m = sqrt(1.0/WM[2,2]);
	
	stdkm = 1.0/sqrt(kappaM);
	stdkf = 1.0/sqrt(kappaF);
	
	for(i in 1:20){
		b1f[i] = wt1f*zf1[i];
		b1f[i+20] = wt12f*zf1[i] + wt2f*zf1[20+i];
		b1m[i] = wt1m*zm1[i];
		b1m[i+20] = wt12m*zm1[i] + wt2m*zm1[20+i];
		b3f[i] = stdkf*zf3[i];
		b3m[i] = stdkm*zm3[i];
	}
	
	
	for(i in 1:240){
		eta[i] = beta[1] + beta[2]*fW[i] + beta[3]*mW[i] + beta[4]*WW[i] + beta[5]*fall[i] 
		+ b1f[f1[i]] + b1m[m1[i]]; 
	}
	
	for(i in 241:360){
		
		eta[i] = beta[1] + beta[2]*fW[i] + beta[3]*mW[i] + beta[4]*WW[i] + beta[5]*fall[i] 
		+ b3f[f2[i-240]] + b3m[m2[i-240]];
	}
	

}


model{
	target += wishart_lpdf(WM | 3.0,fPri);
	target += wishart_lpdf(WF | 3.0,fPri);
	
	
	target += normal_lpdf(zf1 | 0.0,1.0);
	target += normal_lpdf(zm1 | 0.0,1.0);
	target += normal_lpdf(zf3 | 0.0,1.0);
	target += normal_lpdf(zm3 | 0.0,1.0);
	
	target += gamma_lpdf(kappaM | 1.0, 0.622);
	target += gamma_lpdf(kappaF | 1.0, 0.622);
	

	target += bernoulli_logit_lpmf(y | eta);

}

generated quantities{
	real tau1f = 1.0/pow(wt1f,2);
	real tau1m = 1.0/pow(wt1m,2);
	real tau2f = detf/WF[1,1];
	real tau2m = detm/WM[1,1];
	real rhof = -WF[1,2]/sqrt(WF[1,1]*WF[2,2]);
	real rhom = -WM[1,2]/sqrt(WM[1,1]*WM[2,2]);
}





