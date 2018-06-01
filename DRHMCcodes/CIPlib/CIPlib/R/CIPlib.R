
.onAttach<-function(libname,pkgname){

packageStartupMessage(paste("CIPlib, version ",packageVersion("CIPlib"),"\n",
"For usage within a (r)Stan model, do the following steps: \n",
"\n",
"1. Put < #include \"CIPlib.stan\" > \n",
"   inside the functions{} block at the start of your Stan model file. \n \n",
"2. Translate the stan model using stanc_builder() with arguments \n",
"   < allow_undefined=TRUE, isystem=CIP_header_path() >. \n \n",
"3. Compile the stan model using stan_model() with arguments \n",
"   < allow_undefined=TRUE, include=CIP_include() >. \n \n",
"4. Run the model using sampling() \n",
sep=""))
}

#################################################################################
# Get path of header implementing the Stan interfaces
# to c++ code
#################################################################################
CIP_include <- function(){
    return(paste0('\n#include "',normalizePath(paste0(find.package("CIPlib")
    ,'/CIPlib.hpp')), '"\n'))
}

CIP_header_path <- function(){
    return(normalizePath(find.package("CIPlib")))
}



CIP_AR1_psi <- function(omega,T){
    .Call('CIP_AR1_psi',omega,T,PACKAGE='CIPlib')
}


