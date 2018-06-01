#include <cmath>

/*
 Utilities mainly for debugging purposes
 
 */

// overload the Stan AD value extractor
inline double value_of(double arg){return arg;}

/* ############################################################
 
 Functions related to the gamma distribution CIP
 
 ################################################################ */

/*
 Approximation to the g-function associated with
 the gamma distribution default CIP
 */

// computational core
template <class T>
T CIP_gamma_g_core(const T &a){
    if(value_of(a)>0.0){
        return -0.1528257924495051e0 + a
        + 0.1528257924495051e0/
        (0.1e1 + (0.7819628323755627e0 +
                  (0.3868075653216423e0 +
                   (0.1340846511972002e0 +
                    (0.3337571885056357e-1 +
                     (0.6120134586887599e-2 +
                      (0.1011187678928435e-2 +
                       0.2624458484189310e-3 * a)
                      * a) * a) * a) * a) * a) * a);
    } else {
        return -0.3061802078252214e0 + 0.7071067811865475e0 * a
        + 0.3061802078252214e0 /
        (-0.1860384351099194e-2 * sqrt(0.10e1 - a) +
         0.1001860384351099e1 +
         (-0.5672278327747562e0 +
          (0.2083252540656104e0 +
           (-0.5108147181528087e-1 +
            (0.8195932333867585e-2 +
             (-0.1148298720449367e-2 +
              (0.1645337465449455e-3 -
               0.2500498519885696e-4 * a)
              * a) * a) * a) * a) * a) * a);
    }
} // end CIP_gamma_g_core

// first derivative of the gamma g-function
template <class T>
T CIP_gamma_g_deriv_core(const T &a){
    if(value_of(a)>0.0){
        T t1 = a * a;
        T t2 = t1 * a;
        T t3 = t1 * t1;
        T t6 = t3 * t1;
        T t8 = t3 * a;
        T t15 = pow(0.2624458484189310e-3 * t3 * t2 + 0.1011187678928435e-2 *
                    t6 + 0.6120134586887599e-2 * t8 + 0.3337571885056357e-1 *
                    t3 + 0.1340846511972002e0 * t2 + 0.3868075653216423e0 * t1 +
                    0.7819628323755627e0 * a + 0.1e1, 0.2e1);
        return 0.1e1 - 0.1528257924495051e0 / t15 *
        (0.1837120938932517e-2 * t6 + 0.6067126073570610e-2 * t8 +
         0.3060067293443800e-1 * t3 + 0.1335028754022543e0 * t2 +
         0.4022539535916006e0 * t1 + 0.7736151306432846e0 * a +
         0.7819628323755627e0);
    } else {
        T t1 = sqrt(0.2e1);
        T t4 = sqrt(0.10e1 - a);
        T t6 = a * a;
        T t7 = t6 * a;
        T t8 = t6 * t6;
        T t11 = t8 * t6;
        T t13 = t8 * a;
        T t20 = pow(0.1001860384351099e1 - 0.1860384351099194e-2 * t4 -
                    0.2500498519885696e-4 * t8 * t7 + 0.1645337465449455e-3 * t11 -
                    0.1148298720449367e-2 * t13 + 0.8195932333867585e-2 * t8 -
                    0.5108147181528087e-1 * t7 + 0.2083252540656104e0 * t6 -
                    0.5672278327747562e0 * a, 0.2e1);
        return 0.7071067811865475e0 - 0.3061802078252214e0 /
        t20 * (0.9301921755495970e-3 / t4 - 0.1750348963919987e-3 * t11 +
               0.9872024792696730e-3 * t13 - 0.5741493602246835e-2 * t8 +
               0.3278372933547034e-1 * t7 - 0.1532444154458426e0 * t6 +
               0.4166505081312208e0 * a - 0.5672278327747562e0);
    }
    
}


/* ###########################################################################################
 
 Functions related to stationary Gaussian AR(1) model
 
 ############################################################################################# */


/*
 computational core functions for the stationary AR(1) model
 */
double CIP_AR1_psi_intu(const double Z, const double dT){
    double aZ = fabs(Z);
    double emz2 = exp(-2.0*aZ);
    double emz4 = exp(-4.0*aZ);
    double dTm1 = dT-1.0;
    double dTm2 = dT-2.0;
    double dTm3 = dT-3.0;
    double invsqrtT = 1.0/sqrt(dT);
    double t1 = sqrt(1.0 + 2.0*dTm2*emz2 + emz4);
    double logarg = (1.0 + dTm2*emz2 + t1)/(dTm1 + sqrt(2*dTm1));
    double atarg = dTm3*(1.0 - emz2)/(sqrt(2.0*dTm3)*t1);
    
    double t2 = sqrt(2.0*dTm1);
    double t3 = dTm2+emz2;
    double atharg = (dTm1/t2 - t3/t1)/(1.0-dTm1*t3/(t2*t1));
    
    double ret = invsqrtT*(2.0*aZ + log(logarg) + std::atanh(atharg)) + sqrt(2.0*dTm3/dT)*atan(atarg);
    if(Z<0.0){ ret=-ret;}
    return ret;
}

inline double CIP_AR1_psi_u(const double a, const double dT){
    return 2.0*sqrt((1.0 + 2.0*(dT-3.0)/pow(exp(a)+exp(-a),2))/dT);
}
// derivatives of u-function (used to calculate further derivatives of psi(omega) )
inline double CIP_AR1_psi_u_deriv(const double a, const double dT, double &du, double &ddu){
    double expa = exp(a);
    double expma = exp(-a);
    double u = 2.0*sqrt((1.0 + 2.0*(dT-3.0)/pow(expa+expma,2))/dT);
    double t1 = dT*pow(expa+expma,3);
    double t2 = (dT-3.0)*(expa-expma);
    du = -8.0*t2/(u*t1);
    ddu = -64.0*pow(t2,2)/(pow(u,3)*pow(t1,2))
    + 24.0*(dT-3.0)*pow(expa-expma,2)/(u*dT*pow(expa+expma,4))
    - 8.0*(dT-3.0)/(u*dT*pow(expa+expma,2));
    return u;
}

int CIP_AR1_psi_core(const double omega, const int T, double &psi, double &dpsi){
    
    double ao = fabs(omega);
    double dT = static_cast<double>(T);
    
    
    
    // find a crude solution using bisection search
    
    double lb = -1.0e-12; // in case omega=0
    double ub = 0.5*sqrt(dT)*ao;
    double fl = CIP_AR1_psi_intu(lb,dT)-ao;
    double fu = CIP_AR1_psi_intu(ub,dT)-ao;
    double mb,fm;
    int iter = 0;
    
    while(fabs(ub-lb)>1.0e-3 && fabs(fl-fu)>1.0e-3 && iter < 100){
        mb = (iter%2) ? (fl*ub - fu*lb)/(fl-fu) : 0.5*(lb+ub);
        fm = CIP_AR1_psi_intu(mb,dT)-ao;
        if(fm<0.0){
            lb = mb;
            fl = fm;
        } else {
            ub = mb;
            fu = fm;
        }
        iter++;
    }
    
    if(iter==100) {
        psi = 0.0;
        dpsi = 0.0;
        return 1;
    }
    
    // now do some Newton iterations in order to refine the solution
    double res,dres;
    psi = 0.5*(ub+lb);
    iter = 0;
    while(iter<10){
        res = CIP_AR1_psi_intu(psi,dT)-ao;
        dres = CIP_AR1_psi_u(psi,dT);
        if(fabs(res)<1.0e-13) break;
        psi = fmin(ub,fmax(lb,psi-res/dres));
        iter++;
    }
    if(iter==10){
        psi = 0.0;
        dpsi = 0.0;
        return 2;
    }
    // derivative of psi
    dpsi = 1.0/dres;
    // if omega<0, psi is negative
    if(omega<0.0) psi = -psi;
    return 0;
}

/*
 Same as psi_core, but also computes 2nd and 3rd derivs of psi, which
 are used in the Laplace approximation to the default omega prior
 */
int CIP_AR1_psi_core(const double omega,
                     const int T,
                     double &psi,
                     double &dpsi,
                     double &ddpsi,
                     double &dddpsi){
    int retflag = CIP_AR1_psi_core(omega,T,psi,dpsi);
    double u,du,ddu;
    u = CIP_AR1_psi_u_deriv(psi, static_cast<double>(T),du,ddu);
    ddpsi = -du/pow(u,3);
    dddpsi = -(ddu*u - 3.0*pow(du,2))/pow(u,5);
    return retflag;
}



double CIP_AR1_omega_defaultPrior_logkern(const double omega,
                                          const double alpha,
                                          const double beta,
                                          const int T,
                                          double &dlkern,
                                          double &ddlkern){
    double psi,dpsi,ddpsi,dddpsi,phi,vpi,phisq;
    
    CIP_AR1_psi_core(omega,T,psi,dpsi,ddpsi,dddpsi);
    phi = tanh(psi);
    vpi = 0.5*(phi+1.0);
    phisq = pow(phi,2);
    double lkern = (alpha-1.0)*log(vpi) + (beta-1.0)*log(1.0-vpi) +
    log(1.0-phisq) + log(dpsi);
    dlkern = 0.5*(alpha-1.0)*(1.0-phisq)*dpsi/vpi
    - 0.5*(beta-1.0)*(1.0-phisq)*dpsi/(1.0-vpi)
    -2.0*phi*dpsi + ddpsi/dpsi;
    ddlkern = 0.5*(alpha-1.0)*ddpsi*(1.0-phisq)/vpi
    -(alpha-1.0)*pow(dpsi,2)*phi*(1.0-phisq)/vpi
    -0.25*(alpha-1.0)*pow(dpsi,2)*pow(1.0-phisq,2)/pow(vpi,2)
    -0.5*(beta-1.0)*ddpsi*(1.0-phisq)/(1.0-vpi)
    +(beta-1.0)*pow(dpsi,2)*phi*(1.0-phisq)/(1.0-vpi)
    -0.25*(beta-1.0)*pow(dpsi,2)*pow(1.0-phisq,2)/pow(1.0-vpi,2)
    -2.0*pow(dpsi,2)*(1.0-phisq) - 2.0*phi*ddpsi
    +dddpsi/dpsi - pow(ddpsi/dpsi,2);
    return lkern;
}


int CIP_AR1_omega_defaultPrior_prec_core(const double alpha,
                                         const double beta,
                                         const int T,
                                         double &prec){
    // find maximal range for omega where the numerics is still working (tanh(18)<1 in double)
    double omegau = CIP_AR1_psi_intu(18.0, static_cast<double>(T));
    double omegal,omegat,dlkernl,ddlkernl,dlkernu,ddlkernu,dlkernt,ddlkernt,lkern;
    // crude bisection search for maximum
    omegal = -omegau;
    lkern = CIP_AR1_omega_defaultPrior_logkern(omegal,alpha,beta,T,dlkernl,ddlkernl);
    lkern = CIP_AR1_omega_defaultPrior_logkern(omegau,alpha,beta,T,dlkernu,ddlkernu);
    if(dlkernl<0.0 || dlkernu>0.0){
        prec = 1.0;
        return 1;
    }
    int iter = 0;
    while(fabs(omegau-omegal) > 1.0e-3 && iter<100){
        iter++;
        omegat = 0.5*(omegau+omegal);
        lkern = CIP_AR1_omega_defaultPrior_logkern(omegat,alpha,beta,T,dlkernt,ddlkernt);
        if(dlkernt>0.0){
            omegal = omegat;
        } else {
            omegau = omegat;
        }
    }
    if(iter>99){
        prec = 1.0;
        return 2;
    }
    // now do some Newton iterations to refine solution
    omegat = 0.5*(omegau+omegal);
    for(int i=0;i<3;i++){
        lkern = CIP_AR1_omega_defaultPrior_logkern(omegat,alpha,beta,T,dlkernt,ddlkernt);
        omegat = fmin(omegau,fmax(omegal,omegat - dlkernt/ddlkernt));
    }
    lkern = CIP_AR1_omega_defaultPrior_logkern(omegat,alpha,beta,T,dlkernt,ddlkernt);
    prec = -ddlkernt;
    
    return 0;
    
}
