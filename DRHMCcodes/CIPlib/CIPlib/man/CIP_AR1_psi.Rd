\name{CIP_AR1_psi}
\alias{CIP_AR1_psi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computed psi function related to AR(1) model}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CIP_AR1_psi(omega_R, T_R)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{omega_R}{
%%     ~~Describe \code{omega_R} here~~
}
  \item{T_R}{
%%     ~~Describe \code{T_R} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
structure(function (omega_R, T_R) 
.Primitive(".Call")(<pointer: 0x1043423f0>, omega_R, T_R), code = "\n// includes from the plugin\n\n#include <Rcpp.h>\n\n\n#ifndef BEGIN_RCPP\n#define BEGIN_RCPP\n#endif\n\n#ifndef END_RCPP\n#define END_RCPP\n#endif\n\nusing namespace Rcpp;\n\n// user includes\n#include \"/Users/torekleppe/Dropbox/work/hierar_hmc/stan_codes/include/CIPlib_core.hpp\"\n\n// declarations\nextern \"C\" {\nSEXP file6d078efd01a( SEXP omega_R, SEXP T_R) ;\n}\n\n// definition\n\nSEXP file6d078efd01a( SEXP omega_R, SEXP T_R ){\nBEGIN_RCPP\n\ndouble omega = as<double>(omega_R);\nint T = as<int>(T_R);\ndouble psi;\ndouble dpsi;\nCIP_AR1_psi_core(omega,T,psi,dpsi);\nreturn wrap(psi);\n\nEND_RCPP\n}\n\n\n", class = structure("CFunc", package = "inline"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
