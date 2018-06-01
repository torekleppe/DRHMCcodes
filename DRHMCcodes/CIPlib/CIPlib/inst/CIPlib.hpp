#include "CIPlib_core.hpp"

/*
 Stan model interfaces related to gamma distribution.
*/

// evaluate g-function
template <typename T0__>
typename boost::math::tools::promote_args<T0__>::type
CIP_gamma_g(const T0__& a, std::ostream* pstream__){return CIP_gamma_g_core<T0__>(a);}

// evaluate first derivative of g-function
template <typename T0__>
typename boost::math::tools::promote_args<T0__>::type
CIP_gamma_g_deriv(const T0__ &a, std::ostream* pstream__){return CIP_gamma_g_deriv_core<T0__>(a);}


/*
 Stan model interfaces related to AR(1) model.
 */
template <typename T0__>
typename boost::math::tools::promote_args<T0__>::type
CIP_AR1_psi(const T0__ &omega,
          const int &T,
          std::ostream* pstream__){
    // do calculations in double numerics
    double domega = value_of(omega);
    double psi,dpsi;
    int ef = CIP_AR1_psi_core(domega,T,psi,dpsi);
    if (ef>0) {
        std::stringstream errmsg;
        errmsg << "Failed to find solution in CIP_AR1_psi_core, exit flag = " << ef;
        throw std::domain_error(errmsg.str());
    }
    // correct AD using Taylor series trick
    return psi + dpsi*(omega-domega);
}

template <typename T0__>
typename boost::math::tools::promote_args<T0__>::type
CIP_AR1_psi_deriv(const T0__ &omega,
            const int &T,
            std::ostream* pstream__){
    // do calculations in double numerics
    double domega = value_of(omega);
    double psi,dpsi,ddpsi,dddpsi;
    int ef = CIP_AR1_psi_core(domega,T,psi,dpsi,ddpsi,dddpsi);
    if (ef>0) {
        std::stringstream errmsg;
        errmsg << "Failed to find solution in CIP_AR1_psi_core, exit flag = " << ef;
        throw std::domain_error(errmsg.str());
    }
    // correct AD using Taylor series trick
    return dpsi + ddpsi*(omega-domega);
}


template <typename T0__, typename T1__>
typename boost::math::tools::promote_args<T0__, T1__>::type
CIP_AR1_omega_defaultPrior_prec(const T0__ &alpha,
                                const T1__ &beta,
                                const int &T,
                                std::ostream* pstream__){
    double dalpha = value_of(alpha);
    double dbeta = value_of(beta);
    double prec;
    int ef = CIP_AR1_omega_defaultPrior_prec_core(dalpha,dbeta,T,prec);
    if (ef>0) {
        std::stringstream errmsg;
        errmsg << "Failed to find solution in CIP_AR1_omega_defaultPrior_prec_core, exit flag = " << ef;
        throw std::domain_error(errmsg.str());
    }
    return prec;
}




















