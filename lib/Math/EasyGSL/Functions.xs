#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
#include "ppport.h"

#include <gsl/gsl_sf.h>
#include <gsl/gsl_mode.h>

void warn_handler (const char * reason, const char * file, int line, int gsl_errno) {
  warn("Error(gsl-internal): '%s/%s' (%s:%d)", gsl_strerror(gsl_errno), reason, file, line);
  return;
}

double r2double(int rv_, gsl_sf_result result_) {
  if(rv_) warn("Error: rv(%d)='%s' err=%f",rv_,gsl_strerror(rv_),result_.err);
  return result_.val;
}

double e10r2double(int rv_, gsl_sf_result_e10 result_) {
  if(rv_) warn("Error: rv(%d)='%s' err=%f e10=%f",rv_,gsl_strerror(rv_),result_.err,result_.e10);
  return result_.val;
}

size_t arr_ref_size(SV * data, char * warnmsg) {
  size_t rv;
  if ((!SvROK(data)) || (SvTYPE(SvRV(data)) != SVt_PVAV) || ((rv = av_len((AV *)SvRV(data))) < 0)) {
    if (warnmsg) warn(warnmsg);
    return -1;
  }
  return rv+1;
}

MODULE = Math::EasyGSL::Functions	PACKAGE = Math::EasyGSL::Functions

void
GSL_PREC_DOUBLE()
	PPCODE:
		XPUSHs(sv_2mortal(newSViv(GSL_PREC_DOUBLE)));

void
GSL_PREC_SINGLE()
	PPCODE:
		XPUSHs(sv_2mortal(newSViv(GSL_PREC_SINGLE)));

void
GSL_PREC_APPROX()
	PPCODE:
		XPUSHs(sv_2mortal(newSViv(GSL_PREC_APPROX)));

###### generated part - start ######

## GSL function: int gsl_sf_bessel_Yn_array(const int nmin, const int nmax, const double x, double * result_array);
void
sf_bessel_Yn_array(nmin,nmax,x)
		int nmin;
		int nmax;
		double x;
	INIT:
		int rv_;
		double * result_array_;
		int i, result_array_size;
	PPCODE:

		/* allocate memory */
		result_array_size = nmax-nmin+1;
		if (result_array_size<=0) {
                  warn("Warning: sf_bessel_Yn_array() - invalid parameter 'nmax' and/or 'nmin'");
                  XSRETURN_UNDEF;
		}
		result_array_ = malloc(result_array_size * sizeof(double));
		if (!result_array_) {
                  warn("Warning: sf_bessel_Yn_array() - malloc failed");
                  XSRETURN_UNDEF;
		}
		/* do the job */
		rv_ = gsl_sf_bessel_Yn_array(nmin,nmax,x,result_array_);
		if (rv_) warn("Error: rv(%d)",gsl_strerror(rv_));
		for(i=0;i<result_array_size;i++)
		  XPUSHs(sv_2mortal(newSVnv(result_array_[i])));
		/* free memory */
		free(result_array_);

## GSL function: int gsl_sf_bessel_In_array(const int nmin, const int nmax, const double x, double * result_array);
void
sf_bessel_In_array(nmin,nmax,x)
		int nmin;
		int nmax;
		double x;
	INIT:
		int rv_;
		double * result_array_;
		int i, result_array_size;
	PPCODE:

		/* allocate memory */
		result_array_size = nmax-nmin+1;
		if (result_array_size<=0) {
                  warn("Warning: sf_bessel_In_array() - invalid parameter 'nmax' and/or 'nmin'");
                  XSRETURN_UNDEF;
		}
		result_array_ = malloc(result_array_size * sizeof(double));
		if (!result_array_) {
                  warn("Warning: sf_bessel_In_array() - malloc failed");
                  XSRETURN_UNDEF;
		}
		/* do the job */
		rv_ = gsl_sf_bessel_In_array(nmin,nmax,x,result_array_);
		if (rv_) warn("Error: rv(%d)",gsl_strerror(rv_));
		for(i=0;i<result_array_size;i++)
		  XPUSHs(sv_2mortal(newSVnv(result_array_[i])));
		/* free memory */
		free(result_array_);

## GSL function: int gsl_sf_bessel_In_scaled_array(const int nmin, const int nmax, const double x, double * result_array);
void
sf_bessel_In_scaled_array(nmin,nmax,x)
		int nmin;
		int nmax;
		double x;
	INIT:
		int rv_;
		double * result_array_;
		int i, result_array_size;
	PPCODE:

		/* allocate memory */
		result_array_size = nmax-nmin+1;
		if (result_array_size<=0) {
                  warn("Warning: sf_bessel_In_scaled_array() - invalid parameter 'nmax' and/or 'nmin'");
                  XSRETURN_UNDEF;
		}
		result_array_ = malloc(result_array_size * sizeof(double));
		if (!result_array_) {
                  warn("Warning: sf_bessel_In_scaled_array() - malloc failed");
                  XSRETURN_UNDEF;
		}
		/* do the job */
		rv_ = gsl_sf_bessel_In_scaled_array(nmin,nmax,x,result_array_);
		if (rv_) warn("Error: rv(%d)",gsl_strerror(rv_));
		for(i=0;i<result_array_size;i++)
		  XPUSHs(sv_2mortal(newSVnv(result_array_[i])));
		/* free memory */
		free(result_array_);

## GSL function: int gsl_sf_bessel_Kn_array(const int nmin, const int nmax, const double x, double * result_array);
void
sf_bessel_Kn_array(nmin,nmax,x)
		int nmin;
		int nmax;
		double x;
	INIT:
		int rv_;
		double * result_array_;
		int i, result_array_size;
	PPCODE:

		/* allocate memory */
		result_array_size = nmax-nmin+1;
		if (result_array_size<=0) {
                  warn("Warning: sf_bessel_Kn_array() - invalid parameter 'nmax' and/or 'nmin'");
                  XSRETURN_UNDEF;
		}
		result_array_ = malloc(result_array_size * sizeof(double));
		if (!result_array_) {
                  warn("Warning: sf_bessel_Kn_array() - malloc failed");
                  XSRETURN_UNDEF;
		}
		/* do the job */
		rv_ = gsl_sf_bessel_Kn_array(nmin,nmax,x,result_array_);
		if (rv_) warn("Error: rv(%d)",gsl_strerror(rv_));
		for(i=0;i<result_array_size;i++)
		  XPUSHs(sv_2mortal(newSVnv(result_array_[i])));
		/* free memory */
		free(result_array_);

## GSL function: int gsl_sf_bessel_Kn_scaled_array(const int nmin, const int nmax, const double x, double * result_array);
void
sf_bessel_Kn_scaled_array(nmin,nmax,x)
		int nmin;
		int nmax;
		double x;
	INIT:
		int rv_;
		double * result_array_;
		int i, result_array_size;
	PPCODE:

		/* allocate memory */
		result_array_size = nmax-nmin+1;
		if (result_array_size<=0) {
                  warn("Warning: sf_bessel_Kn_scaled_array() - invalid parameter 'nmax' and/or 'nmin'");
                  XSRETURN_UNDEF;
		}
		result_array_ = malloc(result_array_size * sizeof(double));
		if (!result_array_) {
                  warn("Warning: sf_bessel_Kn_scaled_array() - malloc failed");
                  XSRETURN_UNDEF;
		}
		/* do the job */
		rv_ = gsl_sf_bessel_Kn_scaled_array(nmin,nmax,x,result_array_);
		if (rv_) warn("Error: rv(%d)",gsl_strerror(rv_));
		for(i=0;i<result_array_size;i++)
		  XPUSHs(sv_2mortal(newSVnv(result_array_[i])));
		/* free memory */
		free(result_array_);

## GSL function: int gsl_sf_bessel_jl_array(const int lmax, const double x, double * result_array);
void
sf_bessel_jl_array(lmax,x)
		int lmax;
		double x;
	INIT:
		int rv_;
		double * result_array_;
		int i, result_array_size;
	PPCODE:

		/* allocate memory */
		result_array_size = lmax+1;
		if (result_array_size<=0) {
                  warn("Warning: sf_bessel_jl_array() - invalid parameter 'lmax'");
                  XSRETURN_UNDEF;
		}
		result_array_ = malloc(result_array_size * sizeof(double));
		if (!result_array_) {
                  warn("Warning: sf_bessel_jl_array() - malloc failed");
                  XSRETURN_UNDEF;
		}
		/* do the job */
		rv_ = gsl_sf_bessel_jl_array(lmax,x,result_array_);
		if (rv_) warn("Error: rv(%d)",gsl_strerror(rv_));
		for(i=0;i<result_array_size;i++)
		  XPUSHs(sv_2mortal(newSVnv(result_array_[i])));
		/* free memory */
		free(result_array_);

## GSL function: int gsl_sf_bessel_jl_steed_array(const int lmax, const double x, double * jl_x_array);
void
sf_bessel_jl_steed_array(lmax,x)
		int lmax;
		double x;
	INIT:
		int rv_;
		double * result_array_;
		int i, result_array_size;
	PPCODE:

		/* allocate memory */
		result_array_size = lmax+1;
		if (result_array_size<=0) {
                  warn("Warning: sf_bessel_jl_steed_array() - invalid parameter 'lmax'");
                  XSRETURN_UNDEF;
		}
		result_array_ = malloc(result_array_size * sizeof(double));
		if (!result_array_) {
                  warn("Warning: sf_bessel_jl_steed_array() - malloc failed");
                  XSRETURN_UNDEF;
		}
		/* do the job */
		rv_ = gsl_sf_bessel_jl_steed_array(lmax,x,result_array_);
		if (rv_) warn("Error: rv(%d)",gsl_strerror(rv_));
		for(i=0;i<result_array_size;i++)
		  XPUSHs(sv_2mortal(newSVnv(result_array_[i])));
		/* free memory */
		free(result_array_);

## GSL function: int gsl_sf_bessel_yl_array(const int lmax, const double x, double * result_array);
void
sf_bessel_yl_array(lmax,x)
		int lmax;
		double x;
	INIT:
		int rv_;
		double * result_array_;
		int i, result_array_size;
	PPCODE:

		/* allocate memory */
		result_array_size = lmax+1;
		if (result_array_size<=0) {
                  warn("Warning: sf_bessel_yl_array() - invalid parameter 'lmax'");
                  XSRETURN_UNDEF;
		}
		result_array_ = malloc(result_array_size * sizeof(double));
		if (!result_array_) {
                  warn("Warning: sf_bessel_yl_array() - malloc failed");
                  XSRETURN_UNDEF;
		}
		/* do the job */
		rv_ = gsl_sf_bessel_yl_array(lmax,x,result_array_);
		if (rv_) warn("Error: rv(%d)",gsl_strerror(rv_));
		for(i=0;i<result_array_size;i++)
		  XPUSHs(sv_2mortal(newSVnv(result_array_[i])));
		/* free memory */
		free(result_array_);

## GSL function: int gsl_sf_bessel_il_scaled_array(const int lmax, const double x, double * result_array);
void
sf_bessel_il_scaled_array(lmax,x)
		int lmax;
		double x;
	INIT:
		int rv_;
		double * result_array_;
		int i, result_array_size;
	PPCODE:

		/* allocate memory */
		result_array_size = lmax+1;
		if (result_array_size<=0) {
                  warn("Warning: sf_bessel_il_scaled_array() - invalid parameter 'lmax'");
                  XSRETURN_UNDEF;
		}
		result_array_ = malloc(result_array_size * sizeof(double));
		if (!result_array_) {
                  warn("Warning: sf_bessel_il_scaled_array() - malloc failed");
                  XSRETURN_UNDEF;
		}
		/* do the job */
		rv_ = gsl_sf_bessel_il_scaled_array(lmax,x,result_array_);
		if (rv_) warn("Error: rv(%d)",gsl_strerror(rv_));
		for(i=0;i<result_array_size;i++)
		  XPUSHs(sv_2mortal(newSVnv(result_array_[i])));
		/* free memory */
		free(result_array_);

## GSL function: int gsl_sf_bessel_kl_scaled_array(const int lmax, const double x, double * result_array);
void
sf_bessel_kl_scaled_array(lmax,x)
		int lmax;
		double x;
	INIT:
		int rv_;
		double * result_array_;
		int i, result_array_size;
	PPCODE:

		/* allocate memory */
		result_array_size = lmax+1;
		if (result_array_size<=0) {
                  warn("Warning: sf_bessel_kl_scaled_array() - invalid parameter 'lmax'");
                  XSRETURN_UNDEF;
		}
		result_array_ = malloc(result_array_size * sizeof(double));
		if (!result_array_) {
                  warn("Warning: sf_bessel_kl_scaled_array() - malloc failed");
                  XSRETURN_UNDEF;
		}
		/* do the job */
		rv_ = gsl_sf_bessel_kl_scaled_array(lmax,x,result_array_);
		if (rv_) warn("Error: rv(%d)",gsl_strerror(rv_));
		for(i=0;i<result_array_size;i++)
		  XPUSHs(sv_2mortal(newSVnv(result_array_[i])));
		/* free memory */
		free(result_array_);

## GSL function: int gsl_sf_gegenpoly_array(int nmax, double lambda, double x, double * result_array);
void
sf_gegenpoly_array(nmax,lambda,x)
		int nmax;
		double lambda;
		double x;
	INIT:
		int rv_;
		double * result_array_;
		int i, result_array_size;
	PPCODE:

		/* allocate memory */
		result_array_size = nmax+1;
		if (result_array_size<=0) {
                  warn("Warning: sf_gegenpoly_array() - invalid parameter 'nmax'");
                  XSRETURN_UNDEF;
		}
		result_array_ = malloc(result_array_size * sizeof(double));
		if (!result_array_) {
                  warn("Warning: sf_gegenpoly_array() - malloc failed");
                  XSRETURN_UNDEF;
		}
		/* do the job */
		rv_ = gsl_sf_gegenpoly_array(nmax,lambda,x,result_array_);
		if (rv_) warn("Error: rv(%d)",gsl_strerror(rv_));
		for(i=0;i<result_array_size;i++)
		  XPUSHs(sv_2mortal(newSVnv(result_array_[i])));
		/* free memory */
		free(result_array_);

## GSL function: int gsl_sf_legendre_H3d_array(const int lmax, const double lambda, const double eta, double * result_array);
void
sf_legendre_H3d_array(lmax,lambda,eta)
		int lmax;
		double lambda;
		double eta;
	INIT:
		int rv_;
		double * result_array_;
		int i, result_array_size;
	PPCODE:

		/* allocate memory */
		result_array_size = lmax+1;
		if (result_array_size<=0) {
                  warn("Warning: sf_legendre_H3d_array() - invalid parameter 'lmax'");
                  XSRETURN_UNDEF;
		}
		result_array_ = malloc(result_array_size * sizeof(double));
		if (!result_array_) {
                  warn("Warning: sf_legendre_H3d_array() - malloc failed");
                  XSRETURN_UNDEF;
		}
		/* do the job */
		rv_ = gsl_sf_legendre_H3d_array(lmax,lambda,eta,result_array_);
		if (rv_) warn("Error: rv(%d)",gsl_strerror(rv_));
		for(i=0;i<result_array_size;i++)
		  XPUSHs(sv_2mortal(newSVnv(result_array_[i])));
		/* free memory */
		free(result_array_);

## GSL function: int gsl_sf_bessel_Jn_array(int nmin, int nmax, double x, double * result_array);
void
sf_bessel_Jn_array(nmin,nmax,x)
		int nmin;
		int nmax;
		double x;
	INIT:
		int rv_;
		double * result_array_;
		int i, result_array_size;
	PPCODE:

		/* allocate memory */
		result_array_size = nmax-nmin+1;
		if (result_array_size<=0) {
                  warn("Warning: sf_bessel_Jn_array() - invalid parameter 'nmax' and/or 'nmin'");
                  XSRETURN_UNDEF;
		}
		result_array_ = malloc(result_array_size * sizeof(double));
		if (!result_array_) {
                  warn("Warning: sf_bessel_Jn_array() - malloc failed");
                  XSRETURN_UNDEF;
		}
		/* do the job */
		rv_ = gsl_sf_bessel_Jn_array(nmin,nmax,x,result_array_);
		if (rv_) warn("Error: rv(%d)",gsl_strerror(rv_));
		for(i=0;i<result_array_size;i++)
		  XPUSHs(sv_2mortal(newSVnv(result_array_[i])));
		/* free memory */
		free(result_array_);

## GSL function: int gsl_sf_Chi_e(const double x, gsl_sf_result * result);
void
sf_Chi(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_Chi_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_Ci_e(const double x, gsl_sf_result * result);
void
sf_Ci(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_Ci_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_Shi_e(const double x, gsl_sf_result * result);
void
sf_Shi(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_Shi_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_Si_e(const double x, gsl_sf_result * result);
void
sf_Si(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_Si_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_airy_Ai_e(const double x, const gsl_mode_t mode, gsl_sf_result * result);
void
sf_airy_Ai(x,mode)
		double x;
		unsigned int mode;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_airy_Ai_e(x,mode,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_airy_Ai_deriv_e(const double x, gsl_mode_t mode, gsl_sf_result * result);
void
sf_airy_Ai_deriv(x,mode)
		double x;
		unsigned int mode;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_airy_Ai_deriv_e(x,mode,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_airy_Ai_deriv_scaled_e(const double x, gsl_mode_t mode, gsl_sf_result * result);
void
sf_airy_Ai_deriv_scaled(x,mode)
		double x;
		unsigned int mode;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_airy_Ai_deriv_scaled_e(x,mode,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_airy_Ai_scaled_e(const double x, gsl_mode_t mode, gsl_sf_result * result);
void
sf_airy_Ai_scaled(x,mode)
		double x;
		unsigned int mode;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_airy_Ai_scaled_e(x,mode,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_airy_Bi_e(const double x, gsl_mode_t mode, gsl_sf_result * result);
void
sf_airy_Bi(x,mode)
		double x;
		unsigned int mode;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_airy_Bi_e(x,mode,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_airy_Bi_deriv_e(const double x, gsl_mode_t mode, gsl_sf_result * result);
void
sf_airy_Bi_deriv(x,mode)
		double x;
		unsigned int mode;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_airy_Bi_deriv_e(x,mode,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_airy_Bi_deriv_scaled_e(const double x, gsl_mode_t mode, gsl_sf_result * result);
void
sf_airy_Bi_deriv_scaled(x,mode)
		double x;
		unsigned int mode;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_airy_Bi_deriv_scaled_e(x,mode,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_airy_Bi_scaled_e(const double x, gsl_mode_t mode, gsl_sf_result * result);
void
sf_airy_Bi_scaled(x,mode)
		double x;
		unsigned int mode;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_airy_Bi_scaled_e(x,mode,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_airy_zero_Ai_e(unsigned int s, gsl_sf_result * result);
void
sf_airy_zero_Ai(s)
		unsigned int s;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_airy_zero_Ai_e(s,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_airy_zero_Ai_deriv_e(unsigned int s, gsl_sf_result * result);
void
sf_airy_zero_Ai_deriv(s)
		unsigned int s;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_airy_zero_Ai_deriv_e(s,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_airy_zero_Bi_e(unsigned int s, gsl_sf_result * result);
void
sf_airy_zero_Bi(s)
		unsigned int s;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_airy_zero_Bi_e(s,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_airy_zero_Bi_deriv_e(unsigned int s, gsl_sf_result * result);
void
sf_airy_zero_Bi_deriv(s)
		unsigned int s;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_airy_zero_Bi_deriv_e(s,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_angle_restrict_pos_err_e(const double theta, gsl_sf_result * result);
void
sf_angle_restrict_pos_err(theta)
		double theta;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_angle_restrict_pos_err_e(theta,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_angle_restrict_symm_err_e(const double theta, gsl_sf_result * result);
void
sf_angle_restrict_symm_err(theta)
		double theta;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_angle_restrict_symm_err_e(theta,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_atanint_e(const double x, gsl_sf_result * result);
void
sf_atanint(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_atanint_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_bessel_I0_e(const double x, gsl_sf_result * result);
void
sf_bessel_I0(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_bessel_I0_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_bessel_I0_scaled_e(const double x, gsl_sf_result * result);
void
sf_bessel_I0_scaled(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_bessel_I0_scaled_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_bessel_I1_e(const double x, gsl_sf_result * result);
void
sf_bessel_I1(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_bessel_I1_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_bessel_I1_scaled_e(const double x, gsl_sf_result * result);
void
sf_bessel_I1_scaled(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_bessel_I1_scaled_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_bessel_In_e(const int n, const double x, gsl_sf_result * result);
void
sf_bessel_In(n,x)
		int n;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_bessel_In_e(n,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_bessel_In_scaled_e(int n, const double x, gsl_sf_result * result);
void
sf_bessel_In_scaled(n,x)
		int n;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_bessel_In_scaled_e(n,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_bessel_Inu_e(double nu, double x, gsl_sf_result * result);
void
sf_bessel_Inu(nu,x)
		double nu;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_bessel_Inu_e(nu,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_bessel_Inu_scaled_e(double nu, double x, gsl_sf_result * result);
void
sf_bessel_Inu_scaled(nu,x)
		double nu;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_bessel_Inu_scaled_e(nu,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_bessel_J0_e(const double x, gsl_sf_result * result);
void
sf_bessel_J0(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_bessel_J0_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_bessel_J1_e(const double x, gsl_sf_result * result);
void
sf_bessel_J1(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_bessel_J1_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_bessel_Jn_e(int n, double x, gsl_sf_result * result);
void
sf_bessel_Jn(n,x)
		int n;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_bessel_Jn_e(n,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_bessel_Jnu_e(const double nu, const double x, gsl_sf_result * result);
void
sf_bessel_Jnu(nu,x)
		double nu;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_bessel_Jnu_e(nu,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_bessel_K0_e(const double x, gsl_sf_result * result);
void
sf_bessel_K0(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_bessel_K0_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_bessel_K0_scaled_e(const double x, gsl_sf_result * result);
void
sf_bessel_K0_scaled(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_bessel_K0_scaled_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_bessel_K1_e(const double x, gsl_sf_result * result);
void
sf_bessel_K1(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_bessel_K1_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_bessel_K1_scaled_e(const double x, gsl_sf_result * result);
void
sf_bessel_K1_scaled(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_bessel_K1_scaled_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_bessel_Kn_e(const int n, const double x, gsl_sf_result * result);
void
sf_bessel_Kn(n,x)
		int n;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_bessel_Kn_e(n,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_bessel_Kn_scaled_e(int n, const double x, gsl_sf_result * result);
void
sf_bessel_Kn_scaled(n,x)
		int n;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_bessel_Kn_scaled_e(n,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_bessel_Knu_e(const double nu, const double x, gsl_sf_result * result);
void
sf_bessel_Knu(nu,x)
		double nu;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_bessel_Knu_e(nu,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_bessel_Knu_scaled_e(const double nu, const double x, gsl_sf_result * result);
void
sf_bessel_Knu_scaled(nu,x)
		double nu;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_bessel_Knu_scaled_e(nu,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_bessel_Y0_e(const double x, gsl_sf_result * result);
void
sf_bessel_Y0(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_bessel_Y0_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_bessel_Y1_e(const double x, gsl_sf_result * result);
void
sf_bessel_Y1(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_bessel_Y1_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_bessel_Yn_e(int n,const double x, gsl_sf_result * result);
void
sf_bessel_Yn(n,x)
		int n;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_bessel_Yn_e(n,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_bessel_Ynu_e(double nu, double x, gsl_sf_result * result);
void
sf_bessel_Ynu(nu,x)
		double nu;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_bessel_Ynu_e(nu,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_bessel_i0_scaled_e(const double x, gsl_sf_result * result);
void
sf_bessel_i0_scaled(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_bessel_i0_scaled_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_bessel_i1_scaled_e(const double x, gsl_sf_result * result);
void
sf_bessel_i1_scaled(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_bessel_i1_scaled_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_bessel_i2_scaled_e(const double x, gsl_sf_result * result);
void
sf_bessel_i2_scaled(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_bessel_i2_scaled_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_bessel_il_scaled_e(const int l, double x, gsl_sf_result * result);
void
sf_bessel_il_scaled(l,x)
		int l;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_bessel_il_scaled_e(l,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_bessel_j0_e(const double x, gsl_sf_result * result);
void
sf_bessel_j0(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_bessel_j0_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_bessel_j1_e(const double x, gsl_sf_result * result);
void
sf_bessel_j1(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_bessel_j1_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_bessel_j2_e(const double x, gsl_sf_result * result);
void
sf_bessel_j2(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_bessel_j2_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_bessel_jl_e(const int l, const double x, gsl_sf_result * result);
void
sf_bessel_jl(l,x)
		int l;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_bessel_jl_e(l,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_bessel_k0_scaled_e(const double x, gsl_sf_result * result);
void
sf_bessel_k0_scaled(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_bessel_k0_scaled_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_bessel_k1_scaled_e(const double x, gsl_sf_result * result);
void
sf_bessel_k1_scaled(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_bessel_k1_scaled_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_bessel_k2_scaled_e(const double x, gsl_sf_result * result);
void
sf_bessel_k2_scaled(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_bessel_k2_scaled_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_bessel_kl_scaled_e(int l, const double x, gsl_sf_result * result);
void
sf_bessel_kl_scaled(l,x)
		int l;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_bessel_kl_scaled_e(l,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_bessel_lnKnu_e(const double nu, const double x, gsl_sf_result * result);
void
sf_bessel_lnKnu(nu,x)
		double nu;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_bessel_lnKnu_e(nu,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_bessel_y0_e(const double x, gsl_sf_result * result);
void
sf_bessel_y0(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_bessel_y0_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_bessel_y1_e(const double x, gsl_sf_result * result);
void
sf_bessel_y1(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_bessel_y1_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_bessel_y2_e(const double x, gsl_sf_result * result);
void
sf_bessel_y2(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_bessel_y2_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_bessel_yl_e(int l, const double x, gsl_sf_result * result);
void
sf_bessel_yl(l,x)
		int l;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_bessel_yl_e(l,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_bessel_zero_J0_e(unsigned int s, gsl_sf_result * result);
void
sf_bessel_zero_J0(s)
		unsigned int s;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_bessel_zero_J0_e(s,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_bessel_zero_J1_e(unsigned int s, gsl_sf_result * result);
void
sf_bessel_zero_J1(s)
		unsigned int s;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_bessel_zero_J1_e(s,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_bessel_zero_Jnu_e(double nu, unsigned int s, gsl_sf_result * result);
void
sf_bessel_zero_Jnu(nu,s)
		double nu;
		unsigned int s;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_bessel_zero_Jnu_e(nu,s,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_beta_e(const double a, const double b, gsl_sf_result * result);
void
sf_beta(a,b)
		double a;
		double b;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_beta_e(a,b,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_beta_inc_e(const double a, const double b, const double x, gsl_sf_result * result);
void
sf_beta_inc(a,b,x)
		double a;
		double b;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_beta_inc_e(a,b,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_choose_e(unsigned int n, unsigned int m, gsl_sf_result * result);
void
sf_choose(n,m)
		unsigned int n;
		unsigned int m;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_choose_e(n,m,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_clausen_e(double x, gsl_sf_result * result);
void
sf_clausen(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_clausen_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_conicalP_0_e(const double lambda, const double x, gsl_sf_result * result);
void
sf_conicalP_0(lambda,x)
		double lambda;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_conicalP_0_e(lambda,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_conicalP_1_e(const double lambda, const double x, gsl_sf_result * result);
void
sf_conicalP_1(lambda,x)
		double lambda;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_conicalP_1_e(lambda,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_conicalP_cyl_reg_e(const int m, const double lambda, const double x, gsl_sf_result * result);
void
sf_conicalP_cyl_reg(m,lambda,x)
		int m;
		double lambda;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_conicalP_cyl_reg_e(m,lambda,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_conicalP_half_e(const double lambda, const double x, gsl_sf_result * result);
void
sf_conicalP_half(lambda,x)
		double lambda;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_conicalP_half_e(lambda,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_conicalP_mhalf_e(const double lambda, const double x, gsl_sf_result * result);
void
sf_conicalP_mhalf(lambda,x)
		double lambda;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_conicalP_mhalf_e(lambda,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_conicalP_sph_reg_e(const int l, const double lambda, const double x, gsl_sf_result * result);
void
sf_conicalP_sph_reg(l,lambda,x)
		int l;
		double lambda;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_conicalP_sph_reg_e(l,lambda,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_cos_e(double x, gsl_sf_result * result);
void
sf_cos(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_cos_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_cos_err_e(const double x, const double dx, gsl_sf_result * result);
void
sf_cos_err(x,dx)
		double x;
		double dx;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_cos_err_e(x,dx,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_coulomb_CL_e(double L, double eta, gsl_sf_result * result);
void
sf_coulomb_CL(L,eta)
		double L;
		double eta;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_coulomb_CL_e(L,eta,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_coupling_3j_e(int two_ja, int two_jb, int two_jc, int two_ma, int two_mb, int two_mc, gsl_sf_result * result);
void
sf_coupling_3j(two_ja,two_jb,two_jc,two_ma,two_mb,two_mc)
		int two_ja;
		int two_jb;
		int two_jc;
		int two_ma;
		int two_mb;
		int two_mc;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_coupling_3j_e(two_ja,two_jb,two_jc,two_ma,two_mb,two_mc,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_coupling_6j_e(int two_ja, int two_jb, int two_jc, int two_jd, int two_je, int two_jf, gsl_sf_result * result);
void
sf_coupling_6j(two_ja,two_jb,two_jc,two_jd,two_je,two_jf)
		int two_ja;
		int two_jb;
		int two_jc;
		int two_jd;
		int two_je;
		int two_jf;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_coupling_6j_e(two_ja,two_jb,two_jc,two_jd,two_je,two_jf,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_coupling_9j_e(int two_ja, int two_jb, int two_jc, int two_jd, int two_je, int two_jf, int two_jg, int two_jh, int two_ji, gsl_sf_result * result);
void
sf_coupling_9j(two_ja,two_jb,two_jc,two_jd,two_je,two_jf,two_jg,two_jh,two_ji)
		int two_ja;
		int two_jb;
		int two_jc;
		int two_jd;
		int two_je;
		int two_jf;
		int two_jg;
		int two_jh;
		int two_ji;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_coupling_9j_e(two_ja,two_jb,two_jc,two_jd,two_je,two_jf,two_jg,two_jh,two_ji,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_coupling_RacahW_e(int two_ja, int two_jb, int two_jc, int two_jd, int two_je, int two_jf, gsl_sf_result * result);
void
sf_coupling_RacahW(two_ja,two_jb,two_jc,two_jd,two_je,two_jf)
		int two_ja;
		int two_jb;
		int two_jc;
		int two_jd;
		int two_je;
		int two_jf;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_coupling_RacahW_e(two_ja,two_jb,two_jc,two_jd,two_je,two_jf,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_dawson_e(double x, gsl_sf_result * result);
void
sf_dawson(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_dawson_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_debye_1_e(const double x, gsl_sf_result * result);
void
sf_debye_1(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_debye_1_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_debye_2_e(const double x, gsl_sf_result * result);
void
sf_debye_2(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_debye_2_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_debye_3_e(const double x, gsl_sf_result * result);
void
sf_debye_3(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_debye_3_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_debye_4_e(const double x, gsl_sf_result * result);
void
sf_debye_4(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_debye_4_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_debye_5_e(const double x, gsl_sf_result * result);
void
sf_debye_5(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_debye_5_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_debye_6_e(const double x, gsl_sf_result * result);
void
sf_debye_6(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_debye_6_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_dilog_e(const double x, gsl_sf_result * result);
void
sf_dilog(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_dilog_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_doublefact_e(const unsigned int n, gsl_sf_result * result);
void
sf_doublefact(n)
		unsigned int n;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_doublefact_e(n,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_ellint_D_e(double phi, double k, double n, gsl_mode_t mode, gsl_sf_result * result);
void
sf_ellint_D(phi,k,n,mode)
		double phi;
		double k;
		double n;
		unsigned int mode;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_ellint_D_e(phi,k,n,mode,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_ellint_Dcomp_e(double k, gsl_mode_t mode, gsl_sf_result * result);
void
sf_ellint_Dcomp(k,mode)
		double k;
		unsigned int mode;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_ellint_Dcomp_e(k,mode,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_ellint_E_e(double phi, double k, gsl_mode_t mode, gsl_sf_result * result);
void
sf_ellint_E(phi,k,mode)
		double phi;
		double k;
		unsigned int mode;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_ellint_E_e(phi,k,mode,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_ellint_Ecomp_e(double k, gsl_mode_t mode, gsl_sf_result * result);
void
sf_ellint_Ecomp(k,mode)
		double k;
		unsigned int mode;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_ellint_Ecomp_e(k,mode,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_ellint_F_e(double phi, double k, gsl_mode_t mode, gsl_sf_result * result);
void
sf_ellint_F(phi,k,mode)
		double phi;
		double k;
		unsigned int mode;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_ellint_F_e(phi,k,mode,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_ellint_Kcomp_e(double k, gsl_mode_t mode, gsl_sf_result * result);
void
sf_ellint_Kcomp(k,mode)
		double k;
		unsigned int mode;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_ellint_Kcomp_e(k,mode,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_ellint_P_e(double phi, double k, double n, gsl_mode_t mode, gsl_sf_result * result);
void
sf_ellint_P(phi,k,n,mode)
		double phi;
		double k;
		double n;
		unsigned int mode;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_ellint_P_e(phi,k,n,mode,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_ellint_Pcomp_e(double k, double n, gsl_mode_t mode, gsl_sf_result * result);
void
sf_ellint_Pcomp(k,n,mode)
		double k;
		double n;
		unsigned int mode;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_ellint_Pcomp_e(k,n,mode,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_ellint_RC_e(double x, double y, gsl_mode_t mode, gsl_sf_result * result);
void
sf_ellint_RC(x,y,mode)
		double x;
		double y;
		unsigned int mode;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_ellint_RC_e(x,y,mode,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_ellint_RD_e(double x, double y, double z, gsl_mode_t mode, gsl_sf_result * result);
void
sf_ellint_RD(x,y,z,mode)
		double x;
		double y;
		double z;
		unsigned int mode;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_ellint_RD_e(x,y,z,mode,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_ellint_RF_e(double x, double y, double z, gsl_mode_t mode, gsl_sf_result * result);
void
sf_ellint_RF(x,y,z,mode)
		double x;
		double y;
		double z;
		unsigned int mode;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_ellint_RF_e(x,y,z,mode,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_ellint_RJ_e(double x, double y, double z, double p, gsl_mode_t mode, gsl_sf_result * result);
void
sf_ellint_RJ(x,y,z,p,mode)
		double x;
		double y;
		double z;
		double p;
		unsigned int mode;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_ellint_RJ_e(x,y,z,p,mode,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_erf_e(double x, gsl_sf_result * result);
void
sf_erf(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_erf_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_erf_Q_e(double x, gsl_sf_result * result);
void
sf_erf_Q(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_erf_Q_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_erf_Z_e(double x, gsl_sf_result * result);
void
sf_erf_Z(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_erf_Z_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_erfc_e(double x, gsl_sf_result * result);
void
sf_erfc(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_erfc_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_eta_e(const double s, gsl_sf_result * result);
void
sf_eta(s)
		double s;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_eta_e(s,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_eta_int_e(int n, gsl_sf_result * result);
void
sf_eta_int(n)
		int n;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_eta_int_e(n,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_exp_e(const double x, gsl_sf_result * result);
void
sf_exp(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_exp_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_exp_e10_e(const double x, gsl_sf_result_e10 * result);
void
sf_exp_e10(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result_e10 result_;
	PPCODE:
		rv_ = gsl_sf_exp_e10_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(e10r2double(rv_,result_))));

## GSL function: int gsl_sf_exp_err_e(const double x, const double dx, gsl_sf_result * result);
void
sf_exp_err(x,dx)
		double x;
		double dx;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_exp_err_e(x,dx,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_exp_err_e10_e(const double x, const double dx, gsl_sf_result_e10 * result);
void
sf_exp_err_e10(x,dx)
		double x;
		double dx;
	INIT:
		int rv_;
		gsl_sf_result_e10 result_;
	PPCODE:
		rv_ = gsl_sf_exp_err_e10_e(x,dx,&result_);
		XPUSHs(sv_2mortal(newSVnv(e10r2double(rv_,result_))));

## GSL function: int gsl_sf_exp_mult_e(const double x, const double y, gsl_sf_result * result);
void
sf_exp_mult(x,y)
		double x;
		double y;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_exp_mult_e(x,y,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_exp_mult_e10_e(const double x, const double y, gsl_sf_result_e10 * result);
void
sf_exp_mult_e10(x,y)
		double x;
		double y;
	INIT:
		int rv_;
		gsl_sf_result_e10 result_;
	PPCODE:
		rv_ = gsl_sf_exp_mult_e10_e(x,y,&result_);
		XPUSHs(sv_2mortal(newSVnv(e10r2double(rv_,result_))));

## GSL function: int gsl_sf_exp_mult_err_e(const double x, const double dx, const double y, const double dy, gsl_sf_result * result);
void
sf_exp_mult_err(x,dx,y,dy)
		double x;
		double dx;
		double y;
		double dy;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_exp_mult_err_e(x,dx,y,dy,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_exp_mult_err_e10_e(const double x, const double dx, const double y, const double dy, gsl_sf_result_e10 * result);
void
sf_exp_mult_err_e10(x,dx,y,dy)
		double x;
		double dx;
		double y;
		double dy;
	INIT:
		int rv_;
		gsl_sf_result_e10 result_;
	PPCODE:
		rv_ = gsl_sf_exp_mult_err_e10_e(x,dx,y,dy,&result_);
		XPUSHs(sv_2mortal(newSVnv(e10r2double(rv_,result_))));

## GSL function: int gsl_sf_expint_3_e(const double x, gsl_sf_result * result);
void
sf_expint_3(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_expint_3_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_expint_E1_e(const double x, gsl_sf_result * result);
void
sf_expint_E1(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_expint_E1_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_expint_E1_scaled_e(const double x, gsl_sf_result * result);
void
sf_expint_E1_scaled(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_expint_E1_scaled_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_expint_E2_e(const double x, gsl_sf_result * result);
void
sf_expint_E2(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_expint_E2_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_expint_E2_scaled_e(const double x, gsl_sf_result * result);
void
sf_expint_E2_scaled(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_expint_E2_scaled_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_expint_Ei_e(const double x, gsl_sf_result * result);
void
sf_expint_Ei(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_expint_Ei_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_expint_Ei_scaled_e(const double x, gsl_sf_result * result);
void
sf_expint_Ei_scaled(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_expint_Ei_scaled_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_expint_En_e(const int n, const double x, gsl_sf_result * result);
void
sf_expint_En(n,x)
		int n;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_expint_En_e(n,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_expint_En_scaled_e(const int n, const double x, gsl_sf_result * result);
void
sf_expint_En_scaled(n,x)
		int n;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_expint_En_scaled_e(n,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_expm1_e(const double x, gsl_sf_result * result);
void
sf_expm1(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_expm1_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_exprel_e(const double x, gsl_sf_result * result);
void
sf_exprel(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_exprel_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_exprel_2_e(double x, gsl_sf_result * result);
void
sf_exprel_2(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_exprel_2_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_exprel_n_e(const int n, const double x, gsl_sf_result * result);
void
sf_exprel_n(n,x)
		int n;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_exprel_n_e(n,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_exprel_n_CF_e(const double n, const double x, gsl_sf_result * result);
void
sf_exprel_n_CF(n,x)
		double n;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_exprel_n_CF_e(n,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_fact_e(const unsigned int n, gsl_sf_result * result);
void
sf_fact(n)
		unsigned int n;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_fact_e(n,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_fermi_dirac_0_e(const double x, gsl_sf_result * result);
void
sf_fermi_dirac_0(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_fermi_dirac_0_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_fermi_dirac_1_e(const double x, gsl_sf_result * result);
void
sf_fermi_dirac_1(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_fermi_dirac_1_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_fermi_dirac_2_e(const double x, gsl_sf_result * result);
void
sf_fermi_dirac_2(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_fermi_dirac_2_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_fermi_dirac_3half_e(const double x, gsl_sf_result * result);
void
sf_fermi_dirac_3half(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_fermi_dirac_3half_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_fermi_dirac_half_e(const double x, gsl_sf_result * result);
void
sf_fermi_dirac_half(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_fermi_dirac_half_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_fermi_dirac_inc_0_e(const double x, const double b, gsl_sf_result * result);
void
sf_fermi_dirac_inc_0(x,b)
		double x;
		double b;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_fermi_dirac_inc_0_e(x,b,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_fermi_dirac_int_e(const int j, const double x, gsl_sf_result * result);
void
sf_fermi_dirac_int(j,x)
		int j;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_fermi_dirac_int_e(j,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_fermi_dirac_m1_e(const double x, gsl_sf_result * result);
void
sf_fermi_dirac_m1(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_fermi_dirac_m1_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_fermi_dirac_mhalf_e(const double x, gsl_sf_result * result);
void
sf_fermi_dirac_mhalf(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_fermi_dirac_mhalf_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_gamma_e(const double x, gsl_sf_result * result);
void
sf_gamma(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_gamma_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_gamma_inc_e(const double a, const double x, gsl_sf_result * result);
void
sf_gamma_inc(a,x)
		double a;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_gamma_inc_e(a,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_gamma_inc_P_e(const double a, const double x, gsl_sf_result * result);
void
sf_gamma_inc_P(a,x)
		double a;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_gamma_inc_P_e(a,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_gamma_inc_Q_e(const double a, const double x, gsl_sf_result * result);
void
sf_gamma_inc_Q(a,x)
		double a;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_gamma_inc_Q_e(a,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_gammainv_e(const double x, gsl_sf_result * result);
void
sf_gammainv(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_gammainv_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_gammastar_e(const double x, gsl_sf_result * result);
void
sf_gammastar(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_gammastar_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_gegenpoly_1_e(double lambda, double x, gsl_sf_result * result);
void
sf_gegenpoly_1(lambda,x)
		double lambda;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_gegenpoly_1_e(lambda,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_gegenpoly_2_e(double lambda, double x, gsl_sf_result * result);
void
sf_gegenpoly_2(lambda,x)
		double lambda;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_gegenpoly_2_e(lambda,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_gegenpoly_3_e(double lambda, double x, gsl_sf_result * result);
void
sf_gegenpoly_3(lambda,x)
		double lambda;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_gegenpoly_3_e(lambda,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_gegenpoly_n_e(int n, double lambda, double x, gsl_sf_result * result);
void
sf_gegenpoly_n(n,lambda,x)
		int n;
		double lambda;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_gegenpoly_n_e(n,lambda,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_hazard_e(double x, gsl_sf_result * result);
void
sf_hazard(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_hazard_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_hydrogenicR_e(const int n, const int l, const double Z, const double r, gsl_sf_result * result);
void
sf_hydrogenicR(n,l,Z,r)
		int n;
		int l;
		double Z;
		double r;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_hydrogenicR_e(n,l,Z,r,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_hydrogenicR_1_e(const double Z, const double r, gsl_sf_result * result);
void
sf_hydrogenicR_1(Z,r)
		double Z;
		double r;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_hydrogenicR_1_e(Z,r,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_hyperg_0F1_e(double c, double x, gsl_sf_result * result);
void
sf_hyperg_0F1(c,x)
		double c;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_hyperg_0F1_e(c,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_hyperg_1F1_e(const double a, const double b, const double x, gsl_sf_result * result);
void
sf_hyperg_1F1(a,b,x)
		double a;
		double b;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_hyperg_1F1_e(a,b,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_hyperg_1F1_int_e(const int m, const int n, const double x, gsl_sf_result * result);
void
sf_hyperg_1F1_int(m,n,x)
		int m;
		int n;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_hyperg_1F1_int_e(m,n,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_hyperg_2F0_e(const double a, const double b, const double x, gsl_sf_result * result);
void
sf_hyperg_2F0(a,b,x)
		double a;
		double b;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_hyperg_2F0_e(a,b,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_hyperg_2F1_e(double a, double b, const double c, const double x, gsl_sf_result * result);
void
sf_hyperg_2F1(a,b,c,x)
		double a;
		double b;
		double c;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_hyperg_2F1_e(a,b,c,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_hyperg_2F1_conj_e(const double aR, const double aI, const double c, const double x, gsl_sf_result * result);
void
sf_hyperg_2F1_conj(aR,aI,c,x)
		double aR;
		double aI;
		double c;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_hyperg_2F1_conj_e(aR,aI,c,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_hyperg_2F1_conj_renorm_e(const double aR, const double aI, const double c, const double x, gsl_sf_result * result);
void
sf_hyperg_2F1_conj_renorm(aR,aI,c,x)
		double aR;
		double aI;
		double c;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_hyperg_2F1_conj_renorm_e(aR,aI,c,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_hyperg_2F1_renorm_e(const double a, const double b, const double c, const double x, gsl_sf_result * result);
void
sf_hyperg_2F1_renorm(a,b,c,x)
		double a;
		double b;
		double c;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_hyperg_2F1_renorm_e(a,b,c,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_hyperg_U_e(const double a, const double b, const double x, gsl_sf_result * result);
void
sf_hyperg_U(a,b,x)
		double a;
		double b;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_hyperg_U_e(a,b,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_hyperg_U_e10_e(const double a, const double b, const double x, gsl_sf_result_e10 * result);
void
sf_hyperg_U_e10(a,b,x)
		double a;
		double b;
		double x;
	INIT:
		int rv_;
		gsl_sf_result_e10 result_;
	PPCODE:
		rv_ = gsl_sf_hyperg_U_e10_e(a,b,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(e10r2double(rv_,result_))));

## GSL function: int gsl_sf_hyperg_U_int_e(const int m, const int n, const double x, gsl_sf_result * result);
void
sf_hyperg_U_int(m,n,x)
		int m;
		int n;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_hyperg_U_int_e(m,n,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_hyperg_U_int_e10_e(const int m, const int n, const double x, gsl_sf_result_e10 * result);
void
sf_hyperg_U_int_e10(m,n,x)
		int m;
		int n;
		double x;
	INIT:
		int rv_;
		gsl_sf_result_e10 result_;
	PPCODE:
		rv_ = gsl_sf_hyperg_U_int_e10_e(m,n,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(e10r2double(rv_,result_))));

## GSL function: int gsl_sf_hypot_e(const double x, const double y, gsl_sf_result * result);
void
sf_hypot(x,y)
		double x;
		double y;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_hypot_e(x,y,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_hzeta_e(const double s, const double q, gsl_sf_result * result);
void
sf_hzeta(s,q)
		double s;
		double q;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_hzeta_e(s,q,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_laguerre_1_e(const double a, const double x, gsl_sf_result * result);
void
sf_laguerre_1(a,x)
		double a;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_laguerre_1_e(a,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_laguerre_2_e(const double a, const double x, gsl_sf_result * result);
void
sf_laguerre_2(a,x)
		double a;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_laguerre_2_e(a,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_laguerre_3_e(const double a, const double x, gsl_sf_result * result);
void
sf_laguerre_3(a,x)
		double a;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_laguerre_3_e(a,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_laguerre_n_e(const int n, const double a, const double x, gsl_sf_result * result);
void
sf_laguerre_n(n,a,x)
		int n;
		double a;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_laguerre_n_e(n,a,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_lambert_W0_e(double x, gsl_sf_result * result);
void
sf_lambert_W0(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_lambert_W0_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_lambert_Wm1_e(double x, gsl_sf_result * result);
void
sf_lambert_Wm1(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_lambert_Wm1_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_legendre_H3d_e(const int l, const double lambda, const double eta, gsl_sf_result * result);
void
sf_legendre_H3d(l,lambda,eta)
		int l;
		double lambda;
		double eta;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_legendre_H3d_e(l,lambda,eta,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_legendre_H3d_0_e(const double lambda, const double eta, gsl_sf_result * result);
void
sf_legendre_H3d_0(lambda,eta)
		double lambda;
		double eta;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_legendre_H3d_0_e(lambda,eta,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_legendre_H3d_1_e(const double lambda, const double eta, gsl_sf_result * result);
void
sf_legendre_H3d_1(lambda,eta)
		double lambda;
		double eta;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_legendre_H3d_1_e(lambda,eta,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_legendre_P1_e(double x, gsl_sf_result * result);
void
sf_legendre_P1(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_legendre_P1_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_legendre_P2_e(double x, gsl_sf_result * result);
void
sf_legendre_P2(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_legendre_P2_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_legendre_P3_e(double x, gsl_sf_result * result);
void
sf_legendre_P3(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_legendre_P3_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_legendre_Pl_e(const int l, const double x, gsl_sf_result * result);
void
sf_legendre_Pl(l,x)
		int l;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_legendre_Pl_e(l,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_legendre_Plm_e(const int l, const int m, const double x, gsl_sf_result * result);
void
sf_legendre_Plm(l,m,x)
		int l;
		int m;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_legendre_Plm_e(l,m,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_legendre_Q0_e(const double x, gsl_sf_result * result);
void
sf_legendre_Q0(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_legendre_Q0_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_legendre_Q1_e(const double x, gsl_sf_result * result);
void
sf_legendre_Q1(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_legendre_Q1_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_legendre_Ql_e(const int l, const double x, gsl_sf_result * result);
void
sf_legendre_Ql(l,x)
		int l;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_legendre_Ql_e(l,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_legendre_sphPlm_e(const int l, int m, const double x, gsl_sf_result * result);
void
sf_legendre_sphPlm(l,m,x)
		int l;
		int m;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_legendre_sphPlm_e(l,m,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_lnbeta_e(const double a, const double b, gsl_sf_result * result);
void
sf_lnbeta(a,b)
		double a;
		double b;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_lnbeta_e(a,b,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_lnchoose_e(unsigned int n, unsigned int m, gsl_sf_result * result);
void
sf_lnchoose(n,m)
		unsigned int n;
		unsigned int m;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_lnchoose_e(n,m,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_lncosh_e(const double x, gsl_sf_result * result);
void
sf_lncosh(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_lncosh_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_lndoublefact_e(const unsigned int n, gsl_sf_result * result);
void
sf_lndoublefact(n)
		unsigned int n;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_lndoublefact_e(n,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_lnfact_e(const unsigned int n, gsl_sf_result * result);
void
sf_lnfact(n)
		unsigned int n;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_lnfact_e(n,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_lngamma_e(double x, gsl_sf_result * result);
void
sf_lngamma(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_lngamma_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_lnpoch_e(const double a, const double x, gsl_sf_result * result);
void
sf_lnpoch(a,x)
		double a;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_lnpoch_e(a,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_lnsinh_e(const double x, gsl_sf_result * result);
void
sf_lnsinh(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_lnsinh_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_log_e(const double x, gsl_sf_result * result);
void
sf_log(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_log_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_log_1plusx_e(const double x, gsl_sf_result * result);
void
sf_log_1plusx(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_log_1plusx_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_log_1plusx_mx_e(const double x, gsl_sf_result * result);
void
sf_log_1plusx_mx(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_log_1plusx_mx_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_log_abs_e(const double x, gsl_sf_result * result);
void
sf_log_abs(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_log_abs_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_log_erfc_e(double x, gsl_sf_result * result);
void
sf_log_erfc(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_log_erfc_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_multiply_e(const double x, const double y, gsl_sf_result * result);
void
sf_multiply(x,y)
		double x;
		double y;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_multiply_e(x,y,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_multiply_err_e(const double x, const double dx, const double y, const double dy, gsl_sf_result * result);
void
sf_multiply_err(x,dx,y,dy)
		double x;
		double dx;
		double y;
		double dy;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_multiply_err_e(x,dx,y,dy,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_poch_e(const double a, const double x, gsl_sf_result * result);
void
sf_poch(a,x)
		double a;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_poch_e(a,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_pochrel_e(const double a, const double x, gsl_sf_result * result);
void
sf_pochrel(a,x)
		double a;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_pochrel_e(a,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_pow_int_e(double x, int n, gsl_sf_result * result);
void
sf_pow_int(x,n)
		double x;
		int n;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_pow_int_e(x,n,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_psi_e(const double x, gsl_sf_result * result);
void
sf_psi(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_psi_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_psi_1_e(const double x, gsl_sf_result * result);
void
sf_psi_1(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_psi_1_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_psi_1_int_e(const int n, gsl_sf_result * result);
void
sf_psi_1_int(n)
		int n;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_psi_1_int_e(n,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_psi_1piy_e(const double y, gsl_sf_result * result);
void
sf_psi_1piy(y)
		double y;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_psi_1piy_e(y,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_psi_int_e(const int n, gsl_sf_result * result);
void
sf_psi_int(n)
		int n;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_psi_int_e(n,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_psi_n_e(const int n, const double x, gsl_sf_result * result);
void
sf_psi_n(n,x)
		int n;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_psi_n_e(n,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_sin_e(double x, gsl_sf_result * result);
void
sf_sin(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_sin_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_sin_err_e(const double x, const double dx, gsl_sf_result * result);
void
sf_sin_err(x,dx)
		double x;
		double dx;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_sin_err_e(x,dx,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_sinc_e(double x, gsl_sf_result * result);
void
sf_sinc(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_sinc_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_synchrotron_1_e(const double x, gsl_sf_result * result);
void
sf_synchrotron_1(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_synchrotron_1_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_synchrotron_2_e(const double x, gsl_sf_result * result);
void
sf_synchrotron_2(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_synchrotron_2_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_taylorcoeff_e(const int n, const double x, gsl_sf_result * result);
void
sf_taylorcoeff(n,x)
		int n;
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_taylorcoeff_e(n,x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_transport_2_e(const double x, gsl_sf_result * result);
void
sf_transport_2(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_transport_2_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_transport_3_e(const double x, gsl_sf_result * result);
void
sf_transport_3(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_transport_3_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_transport_4_e(const double x, gsl_sf_result * result);
void
sf_transport_4(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_transport_4_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_transport_5_e(const double x, gsl_sf_result * result);
void
sf_transport_5(x)
		double x;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_transport_5_e(x,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_zeta_e(const double s, gsl_sf_result * result);
void
sf_zeta(s)
		double s;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_zeta_e(s,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_zeta_int_e(const int n, gsl_sf_result * result);
void
sf_zeta_int(n)
		int n;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_zeta_int_e(n,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_zetam1_e(const double s, gsl_sf_result * result);
void
sf_zetam1(s)
		double s;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_zetam1_e(s,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

## GSL function: int gsl_sf_zetam1_int_e(const int s, gsl_sf_result * result);
void
sf_zetam1_int(s)
		int s;
	INIT:
		int rv_;
		gsl_sf_result result_;
	PPCODE:
		rv_ = gsl_sf_zetam1_int_e(s,&result_);
		XPUSHs(sv_2mortal(newSVnv(r2double(rv_,result_))));

###### generated part - end ######
