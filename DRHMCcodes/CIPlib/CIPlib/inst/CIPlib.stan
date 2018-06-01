

	/*
		Generic transformed priors
	*/
	real CIP_ExpGammaR_logpdf(real x, real alpha, real beta){
		return gamma_lpdf(exp(x)|alpha,beta) + x;
	}
	

	/*
		Gamma distribution default CIP functions
	*/

	real CIP_gamma_g(real a); // C++ implementation interface
	real CIP_gamma_g_deriv(real a); // C++ implementation interface
	
	real CIP_gamma_a_defaultPrior(real a, real alpha, real beta){
		real g;
		real gderiv;
		g = CIP_gamma_g(a);
		gderiv = CIP_gamma_g_deriv(a);
		return gamma_lpdf(exp(g)|alpha,beta) + g + log(gderiv); 
	}
	
	real CIP_gamma_b_defaultPrior(real b, real alpha, real beta){
		return gamma_lpdf(exp(b)|alpha,beta) + b;
	}
	
	/* 
		Weibull distribution default CIP functions
	*/
	
	// lpdf of Weibull distribution in shape,scale parameterisation
	real CIP_weibull_alpha_lambda_logpdf(real y, real alpha, real lambda){
		return weibull_lpdf(y | alpha,pow(1.526205112,-1.0/alpha)*lambda);
	}
	
	real CIP_weibull_a_defaultPrior(real a, real alpha,real beta){
		return gamma_lpdf(exp(a)|alpha,beta) + a;
	}
	
	real CIP_weibull_b_defaultPrior(real b, real alpha,real beta){
		return gamma_lpdf(exp(b)|alpha,beta) + b;
	}
	
	/*
		stationary Gaussian AR(1) model default CIP functions
	*/
	// psi-function
	real CIP_AR1_psi(real omega, int T); // C++ implementation interface
	real CIP_AR1_psi_deriv(real omega, int T); // C++ implementation interface
	
	
	// approximate precision of the default (transformed Beta) prior for omega
	// Note, only intended for usage in the transformed data section
	real CIP_AR1_omega_defaultPrior_prec(real alpha, real beta, int T); // C++ implementation interface
	
	real CIP_AR1_omega_defaultPrior_logpdf(real omega, real alpha, real beta, int T){
		 real psi = CIP_AR1_psi(omega,T);
		 real dpsi = CIP_AR1_psi_deriv(omega,T);
		 real phi = tanh(psi);
		 real vpi = 0.5*(phi+1.0);
		 return beta_lpdf(vpi | alpha,beta) + log(0.5*(1.0-pow(phi,2))*dpsi);
	}
	
	
	/* ########################################################################
	
	Tri-diagonal Cholesky routines
	
	########################################################################### */
	
	/*
		Cholesky factorization L*L^T of symmetric tridiagonal n x n matrix G, where 
		* diag(G) = [diagFirstLastElem,diagElem,...,diagElem,diagFirstLastElem]
		* first sup/sub-diagonal elements are all equal to offDiagElem
		
		Th routine returns a vector[2*n] L where 
		* L[1:n] is the diagonal of L
		* L[n+1:2*n-1] is the sub-diagonal of L
		* L[2*n] is the log-determinant of L
	*/
	vector CIP_TriDiagChol_const1n(int n, real diagFirstLastElem, real diagElem, real offDiagElem){
		vector[2*n] L;
		real LlogDet;
		// first iteration
		L[1] = sqrt(diagFirstLastElem);
		LlogDet = log(L[1]);
		// iteration 2:n-1
		for ( t in 2:(n-1)){
			L[n+t-1] = offDiagElem/L[t-1];
			L[t] = sqrt(diagElem - pow(L[n+t-1],2));
			LlogDet += log(L[t]);
		}
		// last iteration
		L[2*n-1] = offDiagElem/L[n-1];
		L[n] = sqrt(diagFirstLastElem - pow(L[2*n-1],2));
		LlogDet += log(L[n]);
		// done Cholesky
		
		L[2*n] = LlogDet;
		return(L);
	}
	
	/*
		Cholesky factorization L*L^T of symmetric tridiagonal n x n matrix G, where 
		* diag(G) = diagElem (diagElem has length n)
		* first sup/sub-diagonal elements are all equal to offDiagElem
		
		Th routine returns a vector[2*n] L where 
		* L[1:n] is the diagonal of L
		* L[n+1:2*n-1] is the sub-diagonal of L
		* L[2*n] is the log-determinant of L
	*/
	
	
	vector CIP_TriDiagChol_diag_constod(vector diagElem, real offDiagElem){
		int n = rows(diagElem);
		vector[2*n] L;
		real LlogDet;
		L[1] = sqrt(diagElem[1]);
		LlogDet = log(L[1]);
		for ( t in 2:n){
			L[n+t-1] = offDiagElem/L[t-1];
			L[t] = sqrt(diagElem[t] - pow(L[n+t-1],2));
			LlogDet += log(L[t]);
		}
		L[2*n] = LlogDet;
		return(L);
	}
	
	/*
		Cholesky factorization L*L^T of symmetric tridiagonal n x n matrix G, where 
		* diag(G) = diagElem (diagElem has length n)
		* first sup/sub-diagonal elements are in offDiagElem ( offDiagElem has length n-1)
		
		Th routine returns a vector[2*n] L where 
		* L[1:n] is the diagonal of L
		* L[n+1:2*n-1] is the sub-diagonal of L
		* L[2*n] is the log-determinant of L
	*/
	
	
	vector CIP_TriDiagChol(vector diagElem, vector offDiagElem){
		int n = rows(diagElem);
		vector[2*n] L;
		real LlogDet;
		L[1] = sqrt(diagElem[1]);
		LlogDet = log(L[1]);
		for ( t in 2:n){
			L[n+t-1] = offDiagElem[t-1]/L[t-1];
			L[t] = sqrt(diagElem[t] - pow(L[n+t-1],2));
			LlogDet += log(L[t]);
		}
		L[2*n] = LlogDet;
		return(L);
	}
	
	/*
		Solves L^T x = b for x when L is the output of one of the tridiagonal
		Cholesky factorizations above
	*/
	vector CIP_TriDiagChol_LT_solve(vector L, vector b){
		int n = rows(b);
		vector[n] x;
		// first solve
		x[n] = b[n]/L[n];
		// remaining solves
		for ( tt in 1:(n-1)){
			x[n-tt] = (b[n-tt] - x[n-tt+1]*L[2*n-tt])/L[n-tt];
		}
		return(x);
	}
	
	/*
		Solves L x = b for x when L is the output of one of the tridiagonal
		Cholesky factorizations above
	*/
	vector CIP_TriDiagChol_L_solve(vector L, vector b){
		int n = rows(b);
		vector[n] x;
		// first solve
		x[1] = b[1]/L[1];
		// remaining solves
		for ( i in 2:n){
			x[i] = (b[i] - x[i-1]*L[n+i-1])/L[i];
		}
		return(x);
	}
	
	/*
		Solves L L^T x = G x = b for x when L is the output of one of the tridiagonal
		Cholesky factorizations above
	*/
	vector CIP_TriDiagChol_LLT_solve(vector L, vector b){
		return(CIP_TriDiagChol_LT_solve(L,CIP_TriDiagChol_L_solve(L,b)));
	}
	
	
	


	
	
	
	
	