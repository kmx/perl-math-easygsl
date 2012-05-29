#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
#include "ppport.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

gsl_rng* ref2rng(SV* ref) {
  HV* h;
  SV** s;
  if ((h = (HV*)(SvRV(ref))) == NULL) return NULL;
  if ((s = hv_fetchs(h, "!int!rnghandle", 0)) != NULL) return INT2PTR(gsl_rng*,SvIVX(*s));
  return NULL;
}

size_t arr_ref_size(SV * data, char * warnmsg) {
  size_t rv;
  if ((!SvROK(data)) || (SvTYPE(SvRV(data)) != SVt_PVAV) || ((rv = av_len((AV *)SvRV(data))) < 0)) {
    if (warnmsg) warn(warnmsg);
    return -1;
  }
  return rv+1;
}

MODULE = Math::EasyGSL::Random	PACKAGE = Math::EasyGSL::Random

gsl_rng *
_create(name, env_setup)
        char * name;
        int env_setup;
    CODE:
        const gsl_rng_type ** list;
        const gsl_rng_type * t;
        gsl_rng * r;
        int i;

        /* honor ENV variables */
        if (env_setup) gsl_rng_env_setup();

        /* get rnd_type */
        t = NULL;
        if (name == NULL)
          t = gsl_rng_default;
        else {
          /* get rnd_type by name */
          list = gsl_rng_types_setup();
          for(i=0; list[i] != NULL && t == NULL; i++) {
            if (strcmp(name, list[i]->name) == 0)
              t = list[i];
          }
        }
        if (t == NULL) {
          t = gsl_rng_default;
          warn("Warning: rng type '%s' not found, falling back to default '%s'\n", name, t->name);
        }
        /* alloc rng */
        RETVAL = gsl_rng_alloc(t);

    OUTPUT:
        RETVAL

void
_destroy(r)
        gsl_rng * r;
    CODE:
        gsl_rng_free(r);

const char *
name(self)
        SV * self;
    CODE:
        RETVAL = gsl_rng_name(ref2rng(self));
    OUTPUT:
        RETVAL

void
seed(self,seed)
        SV * self;
        unsigned long seed;
    CODE:
        gsl_rng_set(ref2rng(self), seed);

unsigned long
max(self)
        SV * self;
    CODE:
        RETVAL = gsl_rng_max(ref2rng(self));
    OUTPUT:
        RETVAL

unsigned long
min(self)
        SV * self;
    CODE:
        RETVAL = gsl_rng_min(ref2rng(self));
    OUTPUT:
        RETVAL

void
get_raw(self)
        SV * self;
    PPCODE:
        XPUSHs(sv_2mortal(newSViv(gsl_rng_get(ref2rng(self)))));

void
get_pos(self)
        SV * self;
    PPCODE:
        XPUSHs(sv_2mortal(newSVnv(gsl_rng_uniform_pos(ref2rng(self)))));

void
get(self,...)
        SV * self;
    INIT:
        gsl_rng * r;
    PPCODE:
        r = ref2rng(self);
        if (items==1)
          XPUSHs(sv_2mortal(newSVnv(gsl_rng_uniform(r))));
        else
          XPUSHs(sv_2mortal(newSViv(gsl_rng_uniform_int(r, SvIV(ST(1))))));

void
state(self,...)
        SV * self;
    INIT:
        void * state;
        void * newstate;
        size_t n;
        gsl_rng * r;
    PPCODE:
        r = ref2rng(self);
        state = gsl_rng_state(r);
        n = gsl_rng_size(r);
        if (items==1)
          XPUSHs(sv_2mortal(newSVpv(state,n)));
        else {
          STRLEN len;
          newstate = SvPV(ST(1), len);
          if (len == n)
            memcpy(state, newstate, len);
          else
            warn("Error: state has wrong size %d, expected %d\n", len, n);
        }

###### the following funcs are too complicated to be generated ######

## GSL function: void gsl_ran_sample(const gsl_rng * r, void * dest, size_t k, void * src, size_t n, size_t size);
void
sample(self,k,src)
        SV * self;
        size_t k;
        SV * src;
    INIT:
        size_t * dest;
        size_t * indexes;
        size_t src_size, n, i, size;
        AV * av;
        AV* rv_array;
    PPCODE:
        src_size = arr_ref_size(src, "Warning: sample() - invalid 'src' argument");
        if (src_size<0) XSRETURN_UNDEF;
        indexes = malloc(src_size * sizeof(size_t));
        dest = malloc(k * sizeof(size_t));
        if (!indexes || !dest) {
          warn("Warning: sample() - malloc failed");
          if (indexes) free(indexes);
          if (dest) free(dest);
          XSRETURN_UNDEF;
        }
        for(i=0; i<src_size; i++) indexes[i] = i;
        gsl_ran_sample(ref2rng(self),dest,k,indexes,src_size,sizeof(size_t));
        /* create array to be returned */
        rv_array = (AV *)sv_2mortal((SV *)newAV()); /* new array */
        av_extend(rv_array,k); /* not needed but faster */
        for(i=0, av=(AV*)SvRV(src); i<k; i++) {
          SV ** v = av_fetch(av,dest[i],0);
          av_push(rv_array,newSVsv(*v));
        }
        XPUSHs(sv_2mortal(newRV((SV*)rv_array)));

## GSL function: int gsl_ran_choose(const gsl_rng * r, void * dest, size_t k, void * src, size_t n, size_t size);
void
choose(self,k,src)
        SV * self;
        size_t k;
        SV * src;
    INIT:
        size_t * dest;
        size_t * indexes;
        size_t src_size, n, i, size;
        AV * av;
        AV* rv_array;
    PPCODE:
        src_size = arr_ref_size(src, "Warning: choose() - invalid 'src' argument");
        if (src_size<0) XSRETURN_UNDEF;
        indexes = malloc(src_size * sizeof(size_t));
        dest = malloc(k * sizeof(size_t));
        if (!indexes || !dest) {
          warn("Warning: choose() - malloc failed");
          if (indexes) free(indexes);
          if (dest) free(dest);
          XSRETURN_UNDEF;
        }
        for(i=0; i<src_size; i++) indexes[i] = i;
        gsl_ran_choose(ref2rng(self),dest,k,indexes,src_size,sizeof(size_t));
        /* create array to be returned */
        rv_array = (AV *)sv_2mortal((SV *)newAV()); /* new array */
        av_extend(rv_array,k); /* not needed but faster */
        for(i=0, av=(AV*)SvRV(src); i<k; i++) {
          SV ** v = av_fetch(av,dest[i],0);
          av_push(rv_array,newSVsv(*v));
        }
        XPUSHs(sv_2mortal(newRV((SV*)rv_array)));

## GSL function: void gsl_ran_shuffle(const gsl_rng * r, void * base, size_t nmembm, size_t size);
void
shuffle(self,base)
        SV * self;
        SV * base;
    INIT:
        size_t base_size, i, size;
        size_t * indexes;
        AV * av;
        AV* rv_array;
    PPCODE:
        base_size = arr_ref_size(base, "Warning: shuffle() - invalid 'base' argument");
        if (base_size<0) XSRETURN_UNDEF;
        indexes = malloc(base_size * sizeof(size_t));
        if (!indexes) {
            warn("Warning: shuffle() - malloc failed");
            XSRETURN_UNDEF;
        }
        for(i=0; i<base_size; i++) indexes[i] = i;
        gsl_ran_shuffle(ref2rng(self),indexes,base_size,sizeof(size_t));
        /* create array to be returned */
        rv_array = (AV *)sv_2mortal((SV *)newAV()); /* new array */
        av_extend(rv_array,base_size); /* not needed but faster */
        for(i=0, av=(AV*)SvRV(base); i<base_size; i++) {
          SV ** v = av_fetch(av,indexes[i],0);
          av_push(rv_array,newSVsv(*v));
        }
        XPUSHs(sv_2mortal(newRV((SV*)rv_array)));


## GSL function: void gsl_ran_discrete_preproc(size_t K, const double * P);
## GSL function: size_t gsl_ran_discrete(const gsl_rng * r, const gsl_ran_discrete_t * g);
## GSL function: void gsl_ran_discrete_free(gsl_ran_discrete_t * g);
void
get_discrete(self,P)
        SV * self;
        SV * P;
    INIT:
        gsl_ran_discrete_t * g;
        double * P_;
        size_t rv;
        size_t K;
        AV * av;
        int i;
    PPCODE:
        K = arr_ref_size(P, "Warning: get_discrete() - invalid 'P' argument");
        if (K<0) XSRETURN_UNDEF;
        P_ = malloc(K * sizeof(double));
        if (!P_) {
          warn("Warning: get_discrete() - malloc failed");
          XSRETURN_UNDEF;
        }
        for(i=0, av=(AV*)SvRV(P); i<K; i++) P_[i] = SvNV(*av_fetch(av,i,0));
        g = gsl_ran_discrete_preproc(K,P_);
        rv = gsl_ran_discrete(ref2rng(self),g);        
        XPUSHs(sv_2mortal(newSViv(rv)));        
        gsl_ran_discrete_free(g);
        free(P_);

## GSL function: void gsl_ran_poisson_array(const gsl_rng * r, size_t n, unsigned int array[], double mu);
void
get_poisson_array(self,n,mu)
        SV * self;
        size_t n;
        double mu;
    INIT:
        int i;
        unsigned int * array_;
        AV* rv_array;
    PPCODE:
        array_ = malloc(n * sizeof(unsigned int));
        if (!array_) {
          warn("Warning: get_poisson_array() - malloc failed");
          XSRETURN_UNDEF;
        }
        /* do the job */
        gsl_ran_poisson_array(ref2rng(self),n,array_,mu);
        /* create array to be returned */
        rv_array = (AV *)sv_2mortal((SV *)newAV()); /* new array */
        av_extend(rv_array,n); /* not needed but faster */
        for(i=0; i<n; i++) av_push(rv_array,newSViv(array_[i]));
        XPUSHs(sv_2mortal(newRV((SV*)rv_array)));
        /* free memory */
        free(array_);

## GSL function: void gsl_ran_multinomial(const gsl_rng * r, const size_t K, const unsigned int N, const double p[], unsigned int n[]);
void
get_multinomial(self,N,p)
        SV * self;
        unsigned int N;
        SV * p;
    INIT:
        int i;
        AV * av;
        double * p_;
        size_t K; /* length of p array */
        unsigned int * n_;
        size_t n_size;
        AV* rv_array;
    PPCODE:
        /* check valid ARRAYREF */
        K = arr_ref_size(p, "Warning: get_multinomial() - invalid 'p' argument");
        if (K<0) XSRETURN_UNDEF;
        /* allocate memory */
        p_ = malloc(K * sizeof(double));
        n_ = malloc(K * sizeof(unsigned int));
        if (!p_ || !n_) {
          warn("Warning: get_multinomial() - malloc failed");
          if (p_) free(p_);
          if (n_) free(n_);
          XSRETURN_UNDEF;
        }
        /* copy data */
        for(i=0, av=(AV*)SvRV(p); i<K; i++) p_[i] = SvNV(*av_fetch(av,i,0));	
        /* do the job */
        gsl_ran_multinomial(ref2rng(self),K,N,p_,n_);
        /* create array to be returned */
        rv_array = (AV *)sv_2mortal((SV *)newAV()); /* new array */
        av_extend(rv_array,K); /* not needed but faster */
        for(i=0; i<K; i++) av_push(rv_array,newSViv(n_[i]));
        XPUSHs(sv_2mortal(newRV((SV*)rv_array)));
        /* free memory */
        free(p_);
        free(n_);

## GSL function: void gsl_ran_dirichlet(const gsl_rng * r, const size_t K, const double alpha[], double theta[]);
void
get_dirichlet(self,alpha)
        SV * self;
        SV * alpha;
    INIT:
        int i;
        AV * av;
        double * alpha_;
        size_t alpha_size; /* = K */
        double * theta_;
        AV* rv_array;
    PPCODE:
        /* check valid ARRAYREF */
        alpha_size = arr_ref_size(alpha, "Warning: get_dirichlet() - invalid 'alpha' argument");
        if (alpha_size<0) XSRETURN_UNDEF;
        /* allocate memory */
        alpha_ = malloc(alpha_size * sizeof(double));
        theta_ = malloc(alpha_size * sizeof(double));
        if (!alpha_ || !theta_) {
          warn("Warning: get_dirichlet() - malloc failed");
          if (alpha_) free(alpha_);
          if (theta_) free(theta_);
          XSRETURN_UNDEF;
        }
        /* copy data */
        for(i=0, av=(AV*)SvRV(alpha); i<alpha_size; i++) alpha_[i] = SvNV(*av_fetch(av,i,0));
        /* do the job */
        gsl_ran_dirichlet(ref2rng(self),alpha_size,alpha_,theta_);
        /* create array to be returned */
        rv_array = (AV *)sv_2mortal((SV *)newAV()); /* new array */
        av_extend(rv_array,alpha_size); /* not needed but faster */
        for(i=0; i<alpha_size; i++) av_push(rv_array,newSVnv(theta_[i]));
        XPUSHs(sv_2mortal(newRV((SV*)rv_array)));
        /* free memory */
        free(alpha_);
        free(theta_);

## GSL function: void gsl_ran_bivariate_gaussian(const gsl_rng * r, double sigma_x, double sigma_y, double rho, double *x, double *y);
void
get_bivariate_gaussian(self,sigma_x,sigma_y,rho)
        SV * self;
        double sigma_x;
        double sigma_y;
        double rho;
    INIT:
        double x;
        double y;
    PPCODE:
        gsl_ran_bivariate_gaussian(ref2rng(self),sigma_x,sigma_y,rho,&x,&y);
        XPUSHs(sv_2mortal(newSVnv(x)));
        XPUSHs(sv_2mortal(newSVnv(y)));

## GSL function: void gsl_ran_dir_2d(const gsl_rng * r, double * x, double * y);
void
get_dir_2d(self)
        SV * self;
    INIT:
        double x;
        double y;
    PPCODE:
        gsl_ran_dir_2d(ref2rng(self),&x,&y);
        XPUSHs(sv_2mortal(newSVnv(x)));
        XPUSHs(sv_2mortal(newSVnv(y)));

## GSL function: void gsl_ran_dir_2d_trig_method(const gsl_rng * r, double * x, double * y);
void
get_dir_2d_trig_method(self)
        SV * self;
    INIT:
        double x;
        double y;
    PPCODE:
        gsl_ran_dir_2d_trig_method(ref2rng(self),&x,&y);
        XPUSHs(sv_2mortal(newSVnv(x)));
        XPUSHs(sv_2mortal(newSVnv(y)));

## GSL function: void gsl_ran_dir_3d(const gsl_rng * r, double * x, double * y, double * z);
void
get_dir_3d(self)
        SV * self;
    INIT:
        double x;
        double y;
        double z;
    PPCODE:
        gsl_ran_dir_3d(ref2rng(self),&x,&y,&z);
        XPUSHs(sv_2mortal(newSVnv(x)));
        XPUSHs(sv_2mortal(newSVnv(y)));
        XPUSHs(sv_2mortal(newSVnv(z)));

## GSL function: void gsl_ran_dir_nd(const gsl_rng * r, size_t n, double * x);
void
get_dir_nd(self,n)
        SV * self;
        size_t n;
    INIT:
        double * x;
        size_t i;
        AV* rv_array;
    PPCODE:
        x = malloc(n*sizeof(double));
        if (x !=NULL) {
          gsl_ran_dir_nd(ref2rng(self),n,x);
          for(i=0; i<n; i++) XPUSHs(sv_2mortal(newSVnv(x[i])));
          free(x);
        }

###### generated part - start ######

## GSL function: unsigned int gsl_ran_bernoulli(const gsl_rng * r, double p);
void
get_bernoulli(self,p)
		SV * self;
		double p;
	PPCODE:
		XPUSHs(sv_2mortal(newSViv(gsl_ran_bernoulli(ref2rng(self),p))));

## GSL function: double gsl_ran_beta(const gsl_rng * r, const double a, const double b);
void
get_beta(self,a,b)
		SV * self;
		double a;
		double b;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_beta(ref2rng(self),a,b))));

## GSL function: unsigned int gsl_ran_binomial(const gsl_rng * r, double p, unsigned int n);
void
get_binomial(self,p,n)
		SV * self;
		double p;
		unsigned int n;
	PPCODE:
		XPUSHs(sv_2mortal(newSViv(gsl_ran_binomial(ref2rng(self),p,n))));

## GSL function: unsigned int gsl_ran_binomial_knuth(const gsl_rng * r, double p, unsigned int n);
void
get_binomial_knuth(self,p,n)
		SV * self;
		double p;
		unsigned int n;
	PPCODE:
		XPUSHs(sv_2mortal(newSViv(gsl_ran_binomial_knuth(ref2rng(self),p,n))));

## GSL function: unsigned int gsl_ran_binomial_tpe(const gsl_rng * r, double p, unsigned int n);
void
get_binomial_tpe(self,p,n)
		SV * self;
		double p;
		unsigned int n;
	PPCODE:
		XPUSHs(sv_2mortal(newSViv(gsl_ran_binomial_tpe(ref2rng(self),p,n))));

## GSL function: double gsl_ran_cauchy(const gsl_rng * r, const double a);
void
get_cauchy(self,a)
		SV * self;
		double a;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_cauchy(ref2rng(self),a))));

## GSL function: double gsl_ran_chisq(const gsl_rng * r, const double nu);
void
get_chisq(self,nu)
		SV * self;
		double nu;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_chisq(ref2rng(self),nu))));

## GSL function: double gsl_ran_erlang(const gsl_rng * r, const double a, const double n);
void
get_erlang(self,a,n)
		SV * self;
		double a;
		double n;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_erlang(ref2rng(self),a,n))));

## GSL function: double gsl_ran_exponential(const gsl_rng * r, const double mu);
void
get_exponential(self,mu)
		SV * self;
		double mu;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_exponential(ref2rng(self),mu))));

## GSL function: double gsl_ran_exppow(const gsl_rng * r, const double a, const double b);
void
get_exppow(self,a,b)
		SV * self;
		double a;
		double b;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_exppow(ref2rng(self),a,b))));

## GSL function: double gsl_ran_fdist(const gsl_rng * r, const double nu1, const double nu2);
void
get_fdist(self,nu1,nu2)
		SV * self;
		double nu1;
		double nu2;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_fdist(ref2rng(self),nu1,nu2))));

## GSL function: double gsl_ran_flat(const gsl_rng * r, const double a, const double b);
void
get_flat(self,a,b)
		SV * self;
		double a;
		double b;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_flat(ref2rng(self),a,b))));

## GSL function: double gsl_ran_gamma(const gsl_rng * r, const double a, const double b);
void
get_gamma(self,a,b)
		SV * self;
		double a;
		double b;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_gamma(ref2rng(self),a,b))));

## GSL function: double gsl_ran_gamma_int(const gsl_rng * r, const unsigned int a);
void
get_gamma_int(self,a)
		SV * self;
		unsigned int a;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_gamma_int(ref2rng(self),a))));

## GSL function: double gsl_ran_gamma_knuth(const gsl_rng * r, const double a, const double b);
void
get_gamma_knuth(self,a,b)
		SV * self;
		double a;
		double b;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_gamma_knuth(ref2rng(self),a,b))));

## GSL function: double gsl_ran_gamma_mt(const gsl_rng * r, const double a, const double b);
void
get_gamma_mt(self,a,b)
		SV * self;
		double a;
		double b;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_gamma_mt(ref2rng(self),a,b))));

## GSL function: double gsl_ran_gaussian(const gsl_rng * r, const double sigma);
void
get_gaussian(self,sigma)
		SV * self;
		double sigma;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_gaussian(ref2rng(self),sigma))));

## GSL function: double gsl_ran_gaussian_ratio_method(const gsl_rng * r, const double sigma);
void
get_gaussian_ratio_method(self,sigma)
		SV * self;
		double sigma;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_gaussian_ratio_method(ref2rng(self),sigma))));

## GSL function: double gsl_ran_gaussian_tail(const gsl_rng * r, const double a, const double sigma);
void
get_gaussian_tail(self,a,sigma)
		SV * self;
		double a;
		double sigma;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_gaussian_tail(ref2rng(self),a,sigma))));

## GSL function: double gsl_ran_gaussian_ziggurat(const gsl_rng * r, const double sigma);
void
get_gaussian_ziggurat(self,sigma)
		SV * self;
		double sigma;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_gaussian_ziggurat(ref2rng(self),sigma))));

## GSL function: unsigned int gsl_ran_geometric(const gsl_rng * r, const double p);
void
get_geometric(self,p)
		SV * self;
		double p;
	PPCODE:
		XPUSHs(sv_2mortal(newSViv(gsl_ran_geometric(ref2rng(self),p))));

## GSL function: double gsl_ran_gumbel1(const gsl_rng * r, const double a, const double b);
void
get_gumbel1(self,a,b)
		SV * self;
		double a;
		double b;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_gumbel1(ref2rng(self),a,b))));

## GSL function: double gsl_ran_gumbel2(const gsl_rng * r, const double a, const double b);
void
get_gumbel2(self,a,b)
		SV * self;
		double a;
		double b;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_gumbel2(ref2rng(self),a,b))));

## GSL function: unsigned int gsl_ran_hypergeometric(const gsl_rng * r, unsigned int n1, unsigned int n2, unsigned int t);
void
get_hypergeometric(self,n1,n2,t)
		SV * self;
		unsigned int n1;
		unsigned int n2;
		unsigned int t;
	PPCODE:
		XPUSHs(sv_2mortal(newSViv(gsl_ran_hypergeometric(ref2rng(self),n1,n2,t))));

## GSL function: double gsl_ran_landau(const gsl_rng * r);
void
get_landau(self)
		SV * self;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_landau(ref2rng(self)))));

## GSL function: double gsl_ran_laplace(const gsl_rng * r, const double a);
void
get_laplace(self,a)
		SV * self;
		double a;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_laplace(ref2rng(self),a))));

## GSL function: double gsl_ran_levy(const gsl_rng * r, const double c, const double alpha);
void
get_levy(self,c,alpha)
		SV * self;
		double c;
		double alpha;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_levy(ref2rng(self),c,alpha))));

## GSL function: double gsl_ran_levy_skew(const gsl_rng * r, const double c, const double alpha, const double beta);
void
get_levy_skew(self,c,alpha,beta)
		SV * self;
		double c;
		double alpha;
		double beta;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_levy_skew(ref2rng(self),c,alpha,beta))));

## GSL function: unsigned int gsl_ran_logarithmic(const gsl_rng * r, const double p);
void
get_logarithmic(self,p)
		SV * self;
		double p;
	PPCODE:
		XPUSHs(sv_2mortal(newSViv(gsl_ran_logarithmic(ref2rng(self),p))));

## GSL function: double gsl_ran_logistic(const gsl_rng * r, const double a);
void
get_logistic(self,a)
		SV * self;
		double a;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_logistic(ref2rng(self),a))));

## GSL function: double gsl_ran_lognormal(const gsl_rng * r, const double zeta, const double sigma);
void
get_lognormal(self,zeta,sigma)
		SV * self;
		double zeta;
		double sigma;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_lognormal(ref2rng(self),zeta,sigma))));

## GSL function: unsigned int gsl_ran_negative_binomial(const gsl_rng * r, double p, double n);
void
get_negative_binomial(self,p,n)
		SV * self;
		double p;
		double n;
	PPCODE:
		XPUSHs(sv_2mortal(newSViv(gsl_ran_negative_binomial(ref2rng(self),p,n))));

## GSL function: double gsl_ran_pareto(const gsl_rng * r, double a, const double b);
void
get_pareto(self,a,b)
		SV * self;
		double a;
		double b;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_pareto(ref2rng(self),a,b))));

## GSL function: unsigned int gsl_ran_pascal(const gsl_rng * r, double p, unsigned int n);
void
get_pascal(self,p,n)
		SV * self;
		double p;
		unsigned int n;
	PPCODE:
		XPUSHs(sv_2mortal(newSViv(gsl_ran_pascal(ref2rng(self),p,n))));

## GSL function: unsigned int gsl_ran_poisson(const gsl_rng * r, double mu);
void
get_poisson(self,mu)
		SV * self;
		double mu;
	PPCODE:
		XPUSHs(sv_2mortal(newSViv(gsl_ran_poisson(ref2rng(self),mu))));

## GSL function: double gsl_ran_rayleigh(const gsl_rng * r, const double sigma);
void
get_rayleigh(self,sigma)
		SV * self;
		double sigma;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_rayleigh(ref2rng(self),sigma))));

## GSL function: double gsl_ran_rayleigh_tail(const gsl_rng * r, const double a, const double sigma);
void
get_rayleigh_tail(self,a,sigma)
		SV * self;
		double a;
		double sigma;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_rayleigh_tail(ref2rng(self),a,sigma))));

## GSL function: double gsl_ran_tdist(const gsl_rng * r, const double nu);
void
get_tdist(self,nu)
		SV * self;
		double nu;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_tdist(ref2rng(self),nu))));

## GSL function: double gsl_ran_ugaussian(const gsl_rng * r);
void
get_ugaussian(self)
		SV * self;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_ugaussian(ref2rng(self)))));

## GSL function: double gsl_ran_ugaussian_ratio_method(const gsl_rng * r);
void
get_ugaussian_ratio_method(self)
		SV * self;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_ugaussian_ratio_method(ref2rng(self)))));

## GSL function: double gsl_ran_ugaussian_tail(const gsl_rng * r, const double a);
void
get_ugaussian_tail(self,a)
		SV * self;
		double a;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_ugaussian_tail(ref2rng(self),a))));

## GSL function: double gsl_ran_weibull(const gsl_rng * r, const double a, const double b);
void
get_weibull(self,a,b)
		SV * self;
		double a;
		double b;
	PPCODE:
		XPUSHs(sv_2mortal(newSVnv(gsl_ran_weibull(ref2rng(self),a,b))));

###### generated part - end ######
