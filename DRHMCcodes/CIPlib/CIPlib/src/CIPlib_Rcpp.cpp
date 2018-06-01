#include <Rcpp.h> 
#include "../inst/CIPlib_core.hpp"
using namespace Rcpp;


RcppExport SEXP CIP_AR1_psi( SEXP omega, SEXP T ){
    NumericVector omega_(omega);
    NumericVector T_(T);
    
    int n = omega_.size();
    NumericVector psi_(n);
    NumericVector dpsi_(n);
    
    double psi,dpsi;
    for(int i=0; i<n; i++){
        CIP_AR1_psi_core(omega_[i],static_cast<int>(T_[0]),psi,dpsi);
        psi_[i] = psi;
        dpsi_[i] = dpsi;
    }
    List ret;
    ret["psi"] = psi_;
    ret["dpsi"] = dpsi_;
    ret["omega"] = omega_;
    ret["T"] = T_[0];
    return(ret);
}


