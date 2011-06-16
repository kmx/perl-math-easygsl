#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
#include "ppport.h"

#include <gsl/gsl_randist.h>

void warn_handler (const char * reason, const char * file, int line, int gsl_errno) {
  warn("Error(gsl-internal): '%s/%s' (%s:%d)", gsl_strerror(gsl_errno), reason, file, line);
  return;
}

void croak_handler (const char * reason, const char * file, int line, int gsl_errno) {
  croak("Error(gsl-internal): '%s/%s' (%s:%d)", gsl_strerror(gsl_errno), reason, file, line);
  return;
}

void silent_handler (const char * reason, const char * file, int line, int gsl_errno) {
  return;
}

size_t arr_ref_size(SV * data, char * warnmsg) {
  size_t rv;
  if ((!SvROK(data)) || (SvTYPE(SvRV(data)) != SVt_PVAV) || ((rv = av_len((AV *)SvRV(data))) < 0)) {
    if (warnmsg) warn(warnmsg);
    return -1;
  }
  return rv+1;
}

MODULE = Math::EasyGSL::PDF       PACKAGE = Math::EasyGSL::PDF

BOOT:
/* printf("Hello from the bootstrap - PDF!\n"); / * XXX-FIXME */
gsl_set_error_handler (&warn_handler);

## GSL function: double gsl_ran_dirichlet_pdf(size_t K, const double alpha[], const double theta[]);
void
pdf_dirichlet(K,alpha,theta)
		size_t K;
		SV * alpha;
		SV * theta;
	INIT:
		int i;
		AV * av;
		double * alpha_;
		size_t alpha_size;
		double * theta_;
		size_t theta_size;
	PPCODE:
		/* check valid ARRAYREF */
		alpha_size = arr_ref_size(alpha, "Warning: pdf_dirichlet() - invalid 'alpha' argument");
		theta_size = arr_ref_size(theta, "Warning: pdf_dirichlet() - invalid 'theta' argument");
		if (alpha_size<0 || theta_size<0) XSRETURN_UNDEF;
		/* check size */
		if (alpha_size!=theta_size || alpha_size!=K) {
		  warn("Warning: pdf_dirichlet() - array sizes differ: K(%d), alpha_size(%d), theta_size(%d)", K, alpha_size, theta_size);
		  XSRETURN_UNDEF;
		}
		/* allocate memory */
		alpha_ = malloc(alpha_size * sizeof(double));
		theta_ = malloc(theta_size * sizeof(double));
		if (!alpha_ || !theta_) {
		  warn("Warning: pdf_dirichlet() - malloc failed");
		  if (alpha_) free(alpha_);
		  if (theta_) free(theta_);
		  XSRETURN_UNDEF;
		}
		/* copy data */
		for(i=0, av=(AV*)SvRV(alpha); i<alpha_size; i++) alpha_[i] = SvNV(*av_fetch(av,i,0));	
		for(i=0, av=(AV*)SvRV(theta); i<theta_size; i++) theta_[i] = SvNV(*av_fetch(av,i,0));	
		/* do the job */
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_dirichlet_pdf(K,alpha_,theta_))));
		/* free memory */
		free(alpha_);
		free(theta_);

## GSL function: double gsl_ran_dirichlet_lnpdf(size_t K, const double alpha[], const double theta[]);
void
pdf_dirichlet_ln(K,alpha,theta)
		size_t K;
		SV * alpha;
		SV * theta;
	INIT:
		int i;
		AV * av;
		double * alpha_;
		size_t alpha_size;
		double * theta_;
		size_t theta_size;
	PPCODE:
		/* check valid ARRAYREF */
		alpha_size = arr_ref_size(alpha, "Warning: pdf_dirichlet_ln() - invalid 'alpha' argument");
		theta_size = arr_ref_size(theta, "Warning: pdf_dirichlet_ln() - invalid 'theta' argument");
		if (alpha_size<0 || theta_size<0) XSRETURN_UNDEF;
		/* check size */
		if (alpha_size!=theta_size || alpha_size!=K) {
		  warn("Warning: pdf_dirichlet_ln() - array sizes differ: K(%d), alpha_size(%d), theta_size(%d)", K, alpha_size, theta_size);
		  XSRETURN_UNDEF;
		}
		/* allocate memory */
		alpha_ = malloc(alpha_size * sizeof(double));
		theta_ = malloc(theta_size * sizeof(double));
		if (!alpha_ || !theta_) {
		  warn("Warning: pdf_dirichlet_ln() - malloc failed");
		  if (alpha_) free(alpha_);
		  if (theta_) free(theta_);
		  XSRETURN_UNDEF;
		}
		/* copy data */
		for(i=0, av=(AV*)SvRV(alpha); i<alpha_size; i++) alpha_[i] = SvNV(*av_fetch(av,i,0));	
		for(i=0, av=(AV*)SvRV(theta); i<theta_size; i++) theta_[i] = SvNV(*av_fetch(av,i,0));	
		/* do the job */
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_dirichlet_lnpdf(K,alpha_,theta_))));
		/* free memory */
		free(alpha_);
		free(theta_);

## GSL function: double gsl_ran_multinomial_pdf(size_t K, const double p[], const unsigned int n[]);
void
pdf_multinomial(K,p,n)
		size_t K;
		SV * p;
		SV * n;
	INIT:
		int i;
		AV * av;
		double * p_;
		size_t p_size;
		unsigned int * n_;
		size_t n_size;
	PPCODE:
		/* check valid ARRAYREF */
		p_size = arr_ref_size(p, "Warning: pdf_multinomial() - invalid 'p' argument");
		n_size = arr_ref_size(n, "Warning: pdf_multinomial() - invalid 'n' argument");
		if (p_size<0 || n_size<0) XSRETURN_UNDEF;
		/* check size */
		if (p_size!=n_size || p_size!=K) {
		  warn("Warning: pdf_multinomial() - array sizes differ: K(%d), p_size(%d), n_size(%d)", K, p_size, n_size);
		  XSRETURN_UNDEF;
		}
		/* allocate memory */
		p_ = malloc(p_size * sizeof(double));
		n_ = malloc(n_size * sizeof(unsigned int));
		if (!p_ || !n_) {
		  warn("Warning: pdf_multinomial() - malloc failed");
		  if (p_) free(p_);
		  if (n_) free(n_);
		  XSRETURN_UNDEF;
		}
		/* copy data */
		for(i=0, av=(AV*)SvRV(p); i<p_size; i++) p_[i] = SvNV(*av_fetch(av,i,0));	
		for(i=0, av=(AV*)SvRV(n); i<n_size; i++) n_[i] = SvIV(*av_fetch(av,i,0));	
		/* do the job */
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_multinomial_pdf(K,p_,n_))));
		/* free memory */
		free(p_);
		free(n_);

## GSL function: double gsl_ran_multinomial_lnpdf(size_t K, const double p[], const unsigned int n[]);
void
pdf_multinomial_ln(K,p,n)
		size_t K;
		SV * p;
		SV * n;
	INIT:
		int i;
		AV * av;
		double * p_;
		size_t p_size;
		unsigned int * n_;
		size_t n_size;
	PPCODE:
		/* check valid ARRAYREF */
		p_size = arr_ref_size(p, "Warning: pdf_multinomial_ln() - invalid 'p' argument");
		n_size = arr_ref_size(n, "Warning: pdf_multinomial_ln() - invalid 'n' argument");
		if (p_size<0 || n_size<0) XSRETURN_UNDEF;
		/* check size */
		if (p_size!=n_size || p_size!=K) {
		  warn("Warning: pdf_multinomial_ln() - array sizes differ: K(%d), p_size(%d), n_size(%d)", K, p_size, n_size);
		  XSRETURN_UNDEF;
		}
		/* allocate memory */
		p_ = malloc(p_size * sizeof(double));
		n_ = malloc(n_size * sizeof(unsigned int));
		if (!p_ || !n_) {
		  warn("Warning: pdf_multinomial_ln() - malloc failed");
		  if (p_) free(p_);
		  if (n_) free(n_);
		  XSRETURN_UNDEF;
		}
		/* copy data */
		for(i=0, av=(AV*)SvRV(p); i<p_size; i++) p_[i] = SvNV(*av_fetch(av,i,0));	
		for(i=0, av=(AV*)SvRV(n); i<n_size; i++) n_[i] = SvIV(*av_fetch(av,i,0));	
		/* do the job */
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_multinomial_lnpdf(K,p_,n_))));
		/* free memory */
		free(p_);
		free(n_);


###### generated part - start ######

## GSL function: double gsl_ran_bernoulli_pdf(const unsigned int k, double p);
void
pdf_bernoulli(k,p)
		unsigned int k;
		double p;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_bernoulli_pdf(k,p))));

## GSL function: double gsl_ran_beta_pdf(const double x, const double a, const double b);
void
pdf_beta(x,a,b)
		double x;
		double a;
		double b;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_beta_pdf(x,a,b))));

## GSL function: double gsl_ran_binomial_pdf(const unsigned int k, const double p, const unsigned int n);
void
pdf_binomial(k,p,n)
		unsigned int k;
		double p;
		unsigned int n;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_binomial_pdf(k,p,n))));

## GSL function: double gsl_ran_bivariate_gaussian_pdf(const double x, const double y, const double sigma_x, const double sigma_y, const double rho);
void
pdf_bivariate_gaussian(x,y,sigma_x,sigma_y,rho)
		double x;
		double y;
		double sigma_x;
		double sigma_y;
		double rho;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_bivariate_gaussian_pdf(x,y,sigma_x,sigma_y,rho))));

## GSL function: double gsl_ran_cauchy_pdf(const double x, const double a);
void
pdf_cauchy(x,a)
		double x;
		double a;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_cauchy_pdf(x,a))));

## GSL function: double gsl_ran_chisq_pdf(const double x, const double nu);
void
pdf_chisq(x,nu)
		double x;
		double nu;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_chisq_pdf(x,nu))));

## GSL function: double gsl_ran_erlang_pdf(const double x, const double a, const double n);
void
pdf_erlang(x,a,n)
		double x;
		double a;
		double n;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_erlang_pdf(x,a,n))));

## GSL function: double gsl_ran_exponential_pdf(const double x, const double mu);
void
pdf_exponential(x,mu)
		double x;
		double mu;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_exponential_pdf(x,mu))));

## GSL function: double gsl_ran_exppow_pdf(const double x, const double a, const double b);
void
pdf_exppow(x,a,b)
		double x;
		double a;
		double b;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_exppow_pdf(x,a,b))));

## GSL function: double gsl_ran_fdist_pdf(const double x, const double nu1, const double nu2);
void
pdf_fdist(x,nu1,nu2)
		double x;
		double nu1;
		double nu2;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_fdist_pdf(x,nu1,nu2))));

## GSL function: double gsl_ran_flat_pdf(double x, const double a, const double b);
void
pdf_flat(x,a,b)
		double x;
		double a;
		double b;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_flat_pdf(x,a,b))));

## GSL function: double gsl_ran_gamma_pdf(const double x, const double a, const double b);
void
pdf_gamma(x,a,b)
		double x;
		double a;
		double b;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_gamma_pdf(x,a,b))));

## GSL function: double gsl_ran_gaussian_pdf(const double x, const double sigma);
void
pdf_gaussian(x,sigma)
		double x;
		double sigma;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_gaussian_pdf(x,sigma))));

## GSL function: double gsl_ran_gaussian_tail_pdf(const double x, const double a, const double sigma);
void
pdf_gaussian_tail(x,a,sigma)
		double x;
		double a;
		double sigma;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_gaussian_tail_pdf(x,a,sigma))));

## GSL function: double gsl_ran_geometric_pdf(const unsigned int k, const double p);
void
pdf_geometric(k,p)
		unsigned int k;
		double p;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_geometric_pdf(k,p))));

## GSL function: double gsl_ran_gumbel1_pdf(const double x, const double a, const double b);
void
pdf_gumbel1(x,a,b)
		double x;
		double a;
		double b;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_gumbel1_pdf(x,a,b))));

## GSL function: double gsl_ran_gumbel2_pdf(const double x, const double a, const double b);
void
pdf_gumbel2(x,a,b)
		double x;
		double a;
		double b;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_gumbel2_pdf(x,a,b))));

## GSL function: double gsl_ran_hypergeometric_pdf(const unsigned int k, const unsigned int n1, const unsigned int n2, unsigned int t);
void
pdf_hypergeometric(k,n1,n2,t)
		unsigned int k;
		unsigned int n1;
		unsigned int n2;
		unsigned int t;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_hypergeometric_pdf(k,n1,n2,t))));

## GSL function: double gsl_ran_landau_pdf(const double x);
void
pdf_landau(x)
		double x;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_landau_pdf(x))));

## GSL function: double gsl_ran_laplace_pdf(const double x, const double a);
void
pdf_laplace(x,a)
		double x;
		double a;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_laplace_pdf(x,a))));

## GSL function: double gsl_ran_logarithmic_pdf(const unsigned int k, const double p);
void
pdf_logarithmic(k,p)
		unsigned int k;
		double p;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_logarithmic_pdf(k,p))));

## GSL function: double gsl_ran_logistic_pdf(const double x, const double a);
void
pdf_logistic(x,a)
		double x;
		double a;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_logistic_pdf(x,a))));

## GSL function: double gsl_ran_lognormal_pdf(const double x, const double zeta, const double sigma);
void
pdf_lognormal(x,zeta,sigma)
		double x;
		double zeta;
		double sigma;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_lognormal_pdf(x,zeta,sigma))));

## GSL function: double gsl_ran_negative_binomial_pdf(const unsigned int k, const double p, double n);
void
pdf_negative_binomial(k,p,n)
		unsigned int k;
		double p;
		double n;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_negative_binomial_pdf(k,p,n))));

## GSL function: double gsl_ran_pareto_pdf(const double x, const double a, const double b);
void
pdf_pareto(x,a,b)
		double x;
		double a;
		double b;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_pareto_pdf(x,a,b))));

## GSL function: double gsl_ran_pascal_pdf(const unsigned int k, const double p, unsigned int n);
void
pdf_pascal(k,p,n)
		unsigned int k;
		double p;
		unsigned int n;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_pascal_pdf(k,p,n))));

## GSL function: double gsl_ran_poisson_pdf(const unsigned int k, const double mu);
void
pdf_poisson(k,mu)
		unsigned int k;
		double mu;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_poisson_pdf(k,mu))));

## GSL function: double gsl_ran_rayleigh_pdf(const double x, const double sigma);
void
pdf_rayleigh(x,sigma)
		double x;
		double sigma;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_rayleigh_pdf(x,sigma))));

## GSL function: double gsl_ran_rayleigh_tail_pdf(const double x, const double a, const double sigma);
void
pdf_rayleigh_tail(x,a,sigma)
		double x;
		double a;
		double sigma;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_rayleigh_tail_pdf(x,a,sigma))));

## GSL function: double gsl_ran_tdist_pdf(const double x, const double nu);
void
pdf_tdist(x,nu)
		double x;
		double nu;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_tdist_pdf(x,nu))));

## GSL function: double gsl_ran_ugaussian_pdf(const double x);
void
pdf_ugaussian(x)
		double x;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_ugaussian_pdf(x))));

## GSL function: double gsl_ran_ugaussian_tail_pdf(const double x, const double a);
void
pdf_ugaussian_tail(x,a)
		double x;
		double a;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_ugaussian_tail_pdf(x,a))));

## GSL function: double gsl_ran_weibull_pdf(const double x, const double a, const double b);
void
pdf_weibull(x,a,b)
		double x;
		double a;
		double b;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_weibull_pdf(x,a,b))));

###### generated part - end ######
