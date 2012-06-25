#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
#include "ppport.h"

#include <gsl/gsl_cdf.h>

void warn_handler (const char * reason, const char * file, int line, int gsl_errno) {
  warn("Error(gsl-internal): '%s/%s' (%s:%d)", gsl_strerror(gsl_errno), reason, file, line);
  return;
}

MODULE = Math::EasyGSL::CDF       PACKAGE = Math::EasyGSL::CDF

###### generated part - start ######

## GSL function: double gsl_cdf_beta_P(const double x, const double a, const double b);
void
cdf_beta_P(x,a,b)
        double x;
        double a;
        double b;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_beta_P(x,a,b))));

## GSL function: double gsl_cdf_beta_Pinv(const double P, const double a, const double b);
void
cdf_beta_Pinv(P,a,b)
        double P;
        double a;
        double b;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_beta_Pinv(P,a,b))));

## GSL function: double gsl_cdf_beta_Q(const double x, const double a, const double b);
void
cdf_beta_Q(x,a,b)
        double x;
        double a;
        double b;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_beta_Q(x,a,b))));

## GSL function: double gsl_cdf_beta_Qinv(const double Q, const double a, const double b);
void
cdf_beta_Qinv(Q,a,b)
        double Q;
        double a;
        double b;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_beta_Qinv(Q,a,b))));

## GSL function: double gsl_cdf_binomial_P(const unsigned int k, const double p, const unsigned int n);
void
cdf_binomial_P(k,p,n)
        unsigned int k;
        double p;
        unsigned int n;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_binomial_P(k,p,n))));

## GSL function: double gsl_cdf_binomial_Q(const unsigned int k, const double p, const unsigned int n);
void
cdf_binomial_Q(k,p,n)
        unsigned int k;
        double p;
        unsigned int n;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_binomial_Q(k,p,n))));

## GSL function: double gsl_cdf_cauchy_P(const double x, const double a);
void
cdf_cauchy_P(x,a)
        double x;
        double a;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_cauchy_P(x,a))));

## GSL function: double gsl_cdf_cauchy_Pinv(const double P, const double a);
void
cdf_cauchy_Pinv(P,a)
        double P;
        double a;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_cauchy_Pinv(P,a))));

## GSL function: double gsl_cdf_cauchy_Q(const double x, const double a);
void
cdf_cauchy_Q(x,a)
        double x;
        double a;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_cauchy_Q(x,a))));

## GSL function: double gsl_cdf_cauchy_Qinv(const double Q, const double a);
void
cdf_cauchy_Qinv(Q,a)
        double Q;
        double a;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_cauchy_Qinv(Q,a))));

## GSL function: double gsl_cdf_chisq_P(const double x, const double nu);
void
cdf_chisq_P(x,nu)
        double x;
        double nu;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_chisq_P(x,nu))));

## GSL function: double gsl_cdf_chisq_Pinv(const double P, const double nu);
void
cdf_chisq_Pinv(P,nu)
        double P;
        double nu;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_chisq_Pinv(P,nu))));

## GSL function: double gsl_cdf_chisq_Q(const double x, const double nu);
void
cdf_chisq_Q(x,nu)
        double x;
        double nu;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_chisq_Q(x,nu))));

## GSL function: double gsl_cdf_chisq_Qinv(const double Q, const double nu);
void
cdf_chisq_Qinv(Q,nu)
        double Q;
        double nu;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_chisq_Qinv(Q,nu))));

## GSL function: double gsl_cdf_exponential_P(const double x, const double mu);
void
cdf_exponential_P(x,mu)
        double x;
        double mu;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_exponential_P(x,mu))));

## GSL function: double gsl_cdf_exponential_Pinv(const double P, const double mu);
void
cdf_exponential_Pinv(P,mu)
        double P;
        double mu;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_exponential_Pinv(P,mu))));

## GSL function: double gsl_cdf_exponential_Q(const double x, const double mu);
void
cdf_exponential_Q(x,mu)
        double x;
        double mu;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_exponential_Q(x,mu))));

## GSL function: double gsl_cdf_exponential_Qinv(const double Q, const double mu);
void
cdf_exponential_Qinv(Q,mu)
        double Q;
        double mu;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_exponential_Qinv(Q,mu))));

## GSL function: double gsl_cdf_exppow_P(const double x, const double a, const double b);
void
cdf_exppow_P(x,a,b)
        double x;
        double a;
        double b;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_exppow_P(x,a,b))));

## GSL function: double gsl_cdf_exppow_Q(const double x, const double a, const double b);
void
cdf_exppow_Q(x,a,b)
        double x;
        double a;
        double b;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_exppow_Q(x,a,b))));

## GSL function: double gsl_cdf_fdist_P(const double x, const double nu1, const double nu2);
void
cdf_fdist_P(x,nu1,nu2)
        double x;
        double nu1;
        double nu2;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_fdist_P(x,nu1,nu2))));

## GSL function: double gsl_cdf_fdist_Pinv(const double P, const double nu1, const double nu2);
void
cdf_fdist_Pinv(P,nu1,nu2)
        double P;
        double nu1;
        double nu2;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_fdist_Pinv(P,nu1,nu2))));

## GSL function: double gsl_cdf_fdist_Q(const double x, const double nu1, const double nu2);
void
cdf_fdist_Q(x,nu1,nu2)
        double x;
        double nu1;
        double nu2;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_fdist_Q(x,nu1,nu2))));

## GSL function: double gsl_cdf_fdist_Qinv(const double Q, const double nu1, const double nu2);
void
cdf_fdist_Qinv(Q,nu1,nu2)
        double Q;
        double nu1;
        double nu2;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_fdist_Qinv(Q,nu1,nu2))));

## GSL function: double gsl_cdf_flat_P(const double x, const double a, const double b);
void
cdf_flat_P(x,a,b)
        double x;
        double a;
        double b;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_flat_P(x,a,b))));

## GSL function: double gsl_cdf_flat_Pinv(const double P, const double a, const double b);
void
cdf_flat_Pinv(P,a,b)
        double P;
        double a;
        double b;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_flat_Pinv(P,a,b))));

## GSL function: double gsl_cdf_flat_Q(const double x, const double a, const double b);
void
cdf_flat_Q(x,a,b)
        double x;
        double a;
        double b;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_flat_Q(x,a,b))));

## GSL function: double gsl_cdf_flat_Qinv(const double Q, const double a, const double b);
void
cdf_flat_Qinv(Q,a,b)
        double Q;
        double a;
        double b;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_flat_Qinv(Q,a,b))));

## GSL function: double gsl_cdf_gamma_P(const double x, const double a, const double b);
void
cdf_gamma_P(x,a,b)
        double x;
        double a;
        double b;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_gamma_P(x,a,b))));

## GSL function: double gsl_cdf_gamma_Pinv(const double P, const double a, const double b);
void
cdf_gamma_Pinv(P,a,b)
        double P;
        double a;
        double b;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_gamma_Pinv(P,a,b))));

## GSL function: double gsl_cdf_gamma_Q(const double x, const double a, const double b);
void
cdf_gamma_Q(x,a,b)
        double x;
        double a;
        double b;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_gamma_Q(x,a,b))));

## GSL function: double gsl_cdf_gamma_Qinv(const double Q, const double a, const double b);
void
cdf_gamma_Qinv(Q,a,b)
        double Q;
        double a;
        double b;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_gamma_Qinv(Q,a,b))));

## GSL function: double gsl_cdf_gaussian_P(const double x, const double sigma);
void
cdf_gaussian_P(x,sigma)
        double x;
        double sigma;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_gaussian_P(x,sigma))));

## GSL function: double gsl_cdf_gaussian_Pinv(const double P, const double sigma);
void
cdf_gaussian_Pinv(P,sigma)
        double P;
        double sigma;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_gaussian_Pinv(P,sigma))));

## GSL function: double gsl_cdf_gaussian_Q(const double x, const double sigma);
void
cdf_gaussian_Q(x,sigma)
        double x;
        double sigma;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_gaussian_Q(x,sigma))));

## GSL function: double gsl_cdf_gaussian_Qinv(const double Q, const double sigma);
void
cdf_gaussian_Qinv(Q,sigma)
        double Q;
        double sigma;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_gaussian_Qinv(Q,sigma))));

## GSL function: double gsl_cdf_geometric_P(const unsigned int k, const double p);
void
cdf_geometric_P(k,p)
        unsigned int k;
        double p;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_geometric_P(k,p))));

## GSL function: double gsl_cdf_geometric_Q(const unsigned int k, const double p);
void
cdf_geometric_Q(k,p)
        unsigned int k;
        double p;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_geometric_Q(k,p))));

## GSL function: double gsl_cdf_gumbel1_P(const double x, const double a, const double b);
void
cdf_gumbel1_P(x,a,b)
        double x;
        double a;
        double b;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_gumbel1_P(x,a,b))));

## GSL function: double gsl_cdf_gumbel1_Pinv(const double P, const double a, const double b);
void
cdf_gumbel1_Pinv(P,a,b)
        double P;
        double a;
        double b;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_gumbel1_Pinv(P,a,b))));

## GSL function: double gsl_cdf_gumbel1_Q(const double x, const double a, const double b);
void
cdf_gumbel1_Q(x,a,b)
        double x;
        double a;
        double b;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_gumbel1_Q(x,a,b))));

## GSL function: double gsl_cdf_gumbel1_Qinv(const double Q, const double a, const double b);
void
cdf_gumbel1_Qinv(Q,a,b)
        double Q;
        double a;
        double b;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_gumbel1_Qinv(Q,a,b))));

## GSL function: double gsl_cdf_gumbel2_P(const double x, const double a, const double b);
void
cdf_gumbel2_P(x,a,b)
        double x;
        double a;
        double b;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_gumbel2_P(x,a,b))));

## GSL function: double gsl_cdf_gumbel2_Pinv(const double P, const double a, const double b);
void
cdf_gumbel2_Pinv(P,a,b)
        double P;
        double a;
        double b;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_gumbel2_Pinv(P,a,b))));

## GSL function: double gsl_cdf_gumbel2_Q(const double x, const double a, const double b);
void
cdf_gumbel2_Q(x,a,b)
        double x;
        double a;
        double b;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_gumbel2_Q(x,a,b))));

## GSL function: double gsl_cdf_gumbel2_Qinv(const double Q, const double a, const double b);
void
cdf_gumbel2_Qinv(Q,a,b)
        double Q;
        double a;
        double b;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_gumbel2_Qinv(Q,a,b))));

## GSL function: double gsl_cdf_hypergeometric_P(const unsigned int k, const unsigned int n1, const unsigned int n2, const unsigned int t);
void
cdf_hypergeometric_P(k,n1,n2,t)
        unsigned int k;
        unsigned int n1;
        unsigned int n2;
        unsigned int t;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_hypergeometric_P(k,n1,n2,t))));

## GSL function: double gsl_cdf_hypergeometric_Q(const unsigned int k, const unsigned int n1, const unsigned int n2, const unsigned int t);
void
cdf_hypergeometric_Q(k,n1,n2,t)
        unsigned int k;
        unsigned int n1;
        unsigned int n2;
        unsigned int t;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_hypergeometric_Q(k,n1,n2,t))));

## GSL function: double gsl_cdf_laplace_P(const double x, const double a);
void
cdf_laplace_P(x,a)
        double x;
        double a;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_laplace_P(x,a))));

## GSL function: double gsl_cdf_laplace_Pinv(const double P, const double a);
void
cdf_laplace_Pinv(P,a)
        double P;
        double a;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_laplace_Pinv(P,a))));

## GSL function: double gsl_cdf_laplace_Q(const double x, const double a);
void
cdf_laplace_Q(x,a)
        double x;
        double a;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_laplace_Q(x,a))));

## GSL function: double gsl_cdf_laplace_Qinv(const double Q, const double a);
void
cdf_laplace_Qinv(Q,a)
        double Q;
        double a;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_laplace_Qinv(Q,a))));

## GSL function: double gsl_cdf_logistic_P(const double x, const double a);
void
cdf_logistic_P(x,a)
        double x;
        double a;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_logistic_P(x,a))));

## GSL function: double gsl_cdf_logistic_Pinv(const double P, const double a);
void
cdf_logistic_Pinv(P,a)
        double P;
        double a;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_logistic_Pinv(P,a))));

## GSL function: double gsl_cdf_logistic_Q(const double x, const double a);
void
cdf_logistic_Q(x,a)
        double x;
        double a;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_logistic_Q(x,a))));

## GSL function: double gsl_cdf_logistic_Qinv(const double Q, const double a);
void
cdf_logistic_Qinv(Q,a)
        double Q;
        double a;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_logistic_Qinv(Q,a))));

## GSL function: double gsl_cdf_lognormal_P(const double x, const double zeta, const double sigma);
void
cdf_lognormal_P(x,zeta,sigma)
        double x;
        double zeta;
        double sigma;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_lognormal_P(x,zeta,sigma))));

## GSL function: double gsl_cdf_lognormal_Pinv(const double P, const double zeta, const double sigma);
void
cdf_lognormal_Pinv(P,zeta,sigma)
        double P;
        double zeta;
        double sigma;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_lognormal_Pinv(P,zeta,sigma))));

## GSL function: double gsl_cdf_lognormal_Q(const double x, const double zeta, const double sigma);
void
cdf_lognormal_Q(x,zeta,sigma)
        double x;
        double zeta;
        double sigma;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_lognormal_Q(x,zeta,sigma))));

## GSL function: double gsl_cdf_lognormal_Qinv(const double Q, const double zeta, const double sigma);
void
cdf_lognormal_Qinv(Q,zeta,sigma)
        double Q;
        double zeta;
        double sigma;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_lognormal_Qinv(Q,zeta,sigma))));

## GSL function: double gsl_cdf_negative_binomial_P(const unsigned int k, const double p, const double n);
void
cdf_negative_binomial_P(k,p,n)
        unsigned int k;
        double p;
        double n;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_negative_binomial_P(k,p,n))));

## GSL function: double gsl_cdf_negative_binomial_Q(const unsigned int k, const double p, const double n);
void
cdf_negative_binomial_Q(k,p,n)
        unsigned int k;
        double p;
        double n;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_negative_binomial_Q(k,p,n))));

## GSL function: double gsl_cdf_pareto_P(const double x, const double a, const double b);
void
cdf_pareto_P(x,a,b)
        double x;
        double a;
        double b;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_pareto_P(x,a,b))));

## GSL function: double gsl_cdf_pareto_Pinv(const double P, const double a, const double b);
void
cdf_pareto_Pinv(P,a,b)
        double P;
        double a;
        double b;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_pareto_Pinv(P,a,b))));

## GSL function: double gsl_cdf_pareto_Q(const double x, const double a, const double b);
void
cdf_pareto_Q(x,a,b)
        double x;
        double a;
        double b;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_pareto_Q(x,a,b))));

## GSL function: double gsl_cdf_pareto_Qinv(const double Q, const double a, const double b);
void
cdf_pareto_Qinv(Q,a,b)
        double Q;
        double a;
        double b;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_pareto_Qinv(Q,a,b))));

## GSL function: double gsl_cdf_pascal_P(const unsigned int k, const double p, const unsigned int n);
void
cdf_pascal_P(k,p,n)
        unsigned int k;
        double p;
        unsigned int n;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_pascal_P(k,p,n))));

## GSL function: double gsl_cdf_pascal_Q(const unsigned int k, const double p, const unsigned int n);
void
cdf_pascal_Q(k,p,n)
        unsigned int k;
        double p;
        unsigned int n;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_pascal_Q(k,p,n))));

## GSL function: double gsl_cdf_poisson_P(const unsigned int k, const double mu);
void
cdf_poisson_P(k,mu)
        unsigned int k;
        double mu;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_poisson_P(k,mu))));

## GSL function: double gsl_cdf_poisson_Q(const unsigned int k, const double mu);
void
cdf_poisson_Q(k,mu)
        unsigned int k;
        double mu;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_poisson_Q(k,mu))));

## GSL function: double gsl_cdf_rayleigh_P(const double x, const double sigma);
void
cdf_rayleigh_P(x,sigma)
        double x;
        double sigma;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_rayleigh_P(x,sigma))));

## GSL function: double gsl_cdf_rayleigh_Pinv(const double P, const double sigma);
void
cdf_rayleigh_Pinv(P,sigma)
        double P;
        double sigma;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_rayleigh_Pinv(P,sigma))));

## GSL function: double gsl_cdf_rayleigh_Q(const double x, const double sigma);
void
cdf_rayleigh_Q(x,sigma)
        double x;
        double sigma;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_rayleigh_Q(x,sigma))));

## GSL function: double gsl_cdf_rayleigh_Qinv(const double Q, const double sigma);
void
cdf_rayleigh_Qinv(Q,sigma)
        double Q;
        double sigma;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_rayleigh_Qinv(Q,sigma))));

## GSL function: double gsl_cdf_tdist_P(const double x, const double nu);
void
cdf_tdist_P(x,nu)
        double x;
        double nu;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_tdist_P(x,nu))));

## GSL function: double gsl_cdf_tdist_Pinv(const double P, const double nu);
void
cdf_tdist_Pinv(P,nu)
        double P;
        double nu;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_tdist_Pinv(P,nu))));

## GSL function: double gsl_cdf_tdist_Q(const double x, const double nu);
void
cdf_tdist_Q(x,nu)
        double x;
        double nu;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_tdist_Q(x,nu))));

## GSL function: double gsl_cdf_tdist_Qinv(const double Q, const double nu);
void
cdf_tdist_Qinv(Q,nu)
        double Q;
        double nu;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_tdist_Qinv(Q,nu))));

## GSL function: double gsl_cdf_ugaussian_P(const double x);
void
cdf_ugaussian_P(x)
        double x;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_ugaussian_P(x))));

## GSL function: double gsl_cdf_ugaussian_Pinv(const double P);
void
cdf_ugaussian_Pinv(P)
        double P;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_ugaussian_Pinv(P))));

## GSL function: double gsl_cdf_ugaussian_Q(const double x);
void
cdf_ugaussian_Q(x)
        double x;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_ugaussian_Q(x))));

## GSL function: double gsl_cdf_ugaussian_Qinv(const double Q);
void
cdf_ugaussian_Qinv(Q)
        double Q;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_ugaussian_Qinv(Q))));

## GSL function: double gsl_cdf_weibull_P(const double x, const double a, const double b);
void
cdf_weibull_P(x,a,b)
        double x;
        double a;
        double b;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_weibull_P(x,a,b))));

## GSL function: double gsl_cdf_weibull_Pinv(const double P, const double a, const double b);
void
cdf_weibull_Pinv(P,a,b)
        double P;
        double a;
        double b;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_weibull_Pinv(P,a,b))));

## GSL function: double gsl_cdf_weibull_Q(const double x, const double a, const double b);
void
cdf_weibull_Q(x,a,b)
        double x;
        double a;
        double b;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_weibull_Q(x,a,b))));

## GSL function: double gsl_cdf_weibull_Qinv(const double Q, const double a, const double b);
void
cdf_weibull_Qinv(Q,a,b)
        double Q;
        double a;
        double b;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_cdf_weibull_Qinv(Q,a,b))));

###### generated part - end ######
