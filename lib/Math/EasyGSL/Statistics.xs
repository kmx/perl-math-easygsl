#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
#include "ppport.h"

#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>

void warn_handler (const char * reason, const char * file, int line, int gsl_errno) {
  warn("Error(gsl-internal): '%s/%s' (%s:%d)", gsl_strerror(gsl_errno), reason, file, line);
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

MODULE = Math::EasyGSL::Statistics    PACKAGE = Math::EasyGSL::Statistics

## GSL function: double gsl_stats_median_from_sorted_data(const double sorted_data[], const size_t stride, const size_t n);
void
stats_median(data)
        SV * data;
    INIT:
        int i;
        AV * av;
        double * data_;
        size_t data_size;
    PPCODE:
        /* check valid ARRAYREF */
        data_size = arr_ref_size(data, "Warning: stats_median() - invalid 'data' argument");
        if (data_size<0) XSRETURN_UNDEF;
        /* allocate memory */
        data_ = malloc(data_size * sizeof(double));
        if (!data_) {
          warn("Warning: stats_median() - malloc failed");
          XSRETURN_UNDEF;
        }
        /* copy data */
        for(i=0, av=(AV*)SvRV(data); i<data_size; i++) data_[i] = SvNV(*av_fetch(av,i,0));
        /* do the job */
        gsl_sort(data_, 1, data_size);
        XPUSHs(sv_2mortal(newSVnv(gsl_stats_median_from_sorted_data(data_,1,data_size))));
        /* free memory */
        free(data_);

## GSL function: void gsl_stats_minmax(double * min, double * max, const double data[], const size_t stride, const size_t n);
void
stats_minmax(data)
        SV * data;
    INIT:
        int i;
        AV * av;
        double min_;
        double max_;
        double * data_;
        size_t data_size;
    PPCODE:
        /* check valid ARRAYREF */
        data_size = arr_ref_size(data, "Warning: stats_minmax() - invalid 'data' argument");
        if (data_size<0) XSRETURN_UNDEF;
        /* allocate memory */
        data_ = malloc(data_size * sizeof(double));
        if (!data_) {
          warn("Warning: stats_minmax() - malloc failed");
          XSRETURN_UNDEF;
        }
        /* copy data */
        for(i=0, av=(AV*)SvRV(data); i<data_size; i++) data_[i] = SvNV(*av_fetch(av,i,0));
        /* do the job */
        gsl_stats_minmax(&min_,&max_,data_,1,data_size);
        XPUSHs(sv_2mortal(newSVnv(min_)));
        XPUSHs(sv_2mortal(newSVnv(max_)));
        /* free memory */
        free(data_);

## GSL function: void gsl_stats_minmax_index(size_t * min_index, size_t * max_index, const double data[], const size_t stride, const size_t n);
void
stats_minmax_index(data)
        SV * data;
    INIT:
        int i;
        AV * av;
        size_t min_index_;
        size_t max_index_;
        double * data_;
        size_t data_size;
    PPCODE:
        /* check valid ARRAYREF */
        data_size = arr_ref_size(data, "Warning: stats_minmax_index() - invalid 'data' argument");
        if (data_size<0) XSRETURN_UNDEF;
        /* allocate memory */
        data_ = malloc(data_size * sizeof(double));
        if (!data_) {
          warn("Warning: stats_minmax_index() - malloc failed");
          XSRETURN_UNDEF;
        }
        /* copy data */
        for(i=0, av=(AV*)SvRV(data); i<data_size; i++) data_[i] = SvNV(*av_fetch(av,i,0));
        /* do the job */
        gsl_stats_minmax_index(&min_index_,&max_index_,data_,1,data_size);
        XPUSHs(sv_2mortal(newSVnv(min_index_)));
        XPUSHs(sv_2mortal(newSVnv(max_index_)));
        /* free memory */
        free(data_);

## GSL function: double gsl_stats_quantile_from_sorted_data(const double sorted_data[], const size_t stride, const size_t n, const double f);
void
stats_quantile(data,f)
        SV * data;
        double f;
    INIT:
        int i;
        AV * av;
        double * data_;
        size_t data_size;
    PPCODE:
        /* check valid ARRAYREF */
        data_size = arr_ref_size(data, "Warning: stats_quantile() - invalid 'data' argument");
        if (data_size<0) XSRETURN_UNDEF;
        /* allocate memory */
        data_ = malloc(data_size * sizeof(double));
        if (!data_) {
          warn("Warning: stats_quantile() - malloc failed");
          XSRETURN_UNDEF;
        }
        /* copy data */
        for(i=0, av=(AV*)SvRV(data); i<data_size; i++) data_[i] = SvNV(*av_fetch(av,i,0));
        /* do the job */
        XPUSHs(sv_2mortal(newSVnv(gsl_stats_quantile_from_sorted_data(data_,1,data_size,f))));
        /* free memory */
        free(data_);

###### generated part - start ######

## GSL function: double gsl_stats_absdev(const double data[], const size_t stride, const size_t n);
void
stats_absdev(data)
        SV * data;
    INIT:
        int i;
        AV * av;
        double * data_;
        size_t data_size;
    PPCODE:
        /* check valid ARRAYREF */
        data_size = arr_ref_size(data, "Warning: stats_absdev() - invalid 'data' argument");
        if (data_size<0) XSRETURN_UNDEF;
        /* allocate memory */
        data_ = malloc(data_size * sizeof(double));
        if (!data_) {
          warn("Warning: stats_absdev() - malloc failed");
          XSRETURN_UNDEF;
        }
        /* copy data */
        for(i=0, av=(AV*)SvRV(data); i<data_size; i++) data_[i] = SvNV(*av_fetch(av,i,0));
        /* do the job */
        XPUSHs(sv_2mortal(newSVnv(gsl_stats_absdev(data_,1,data_size))));
        /* free memory */
        free(data_);

## GSL function: double gsl_stats_absdev_m(const double data[], const size_t stride, const size_t n, const double mean);
void
stats_absdev_m(data,mean)
        SV * data;
        double mean;
    INIT:
        int i;
        AV * av;
        double * data_;
        size_t data_size;
    PPCODE:
        /* check valid ARRAYREF */
        data_size = arr_ref_size(data, "Warning: stats_absdev_m() - invalid 'data' argument");
        if (data_size<0) XSRETURN_UNDEF;
        /* allocate memory */
        data_ = malloc(data_size * sizeof(double));
        if (!data_) {
          warn("Warning: stats_absdev_m() - malloc failed");
          XSRETURN_UNDEF;
        }
        /* copy data */
        for(i=0, av=(AV*)SvRV(data); i<data_size; i++) data_[i] = SvNV(*av_fetch(av,i,0));
        /* do the job */
        XPUSHs(sv_2mortal(newSVnv(gsl_stats_absdev_m(data_,1,data_size,mean))));
        /* free memory */
        free(data_);

## GSL function: double gsl_stats_correlation(const double data1[], const size_t stride1,const double data2[], const size_t stride2, const size_t n1);
void
stats_correlation(data1,data2)
        SV * data1;
        SV * data2;
    INIT:
        int i;
        AV * av;
        double * data1_;
        size_t data1_size;
        double * data2_;
        size_t data2_size;
    PPCODE:
        /* check valid ARRAYREF */
        data1_size = arr_ref_size(data1, "Warning: stats_correlation() - invalid 'data1' argument");
        data2_size = arr_ref_size(data2, "Warning: stats_correlation() - invalid 'data2' argument");
        if (data1_size<0 || data2_size<0) XSRETURN_UNDEF;
        /* check size */
        if (data1_size!=data2_size) {
          warn("Warning: stats_correlation() - array sizes differ: data1_size(%d), data2_size(%d)", data1_size, data2_size);
          XSRETURN_UNDEF;
        }
        /* allocate memory */
        data1_ = malloc(data1_size * sizeof(double));
        data2_ = malloc(data2_size * sizeof(double));
        if (!data1_ || !data2_) {
          warn("Warning: stats_correlation() - malloc failed");
          if (data1_) free(data1_);
          if (data2_) free(data2_);
          XSRETURN_UNDEF;
        }
        /* copy data */
        for(i=0, av=(AV*)SvRV(data1); i<data1_size; i++) data1_[i] = SvNV(*av_fetch(av,i,0));
        for(i=0, av=(AV*)SvRV(data2); i<data2_size; i++) data2_[i] = SvNV(*av_fetch(av,i,0));
        /* do the job */
        XPUSHs(sv_2mortal(newSVnv(gsl_stats_correlation(data1_,1,data2_,1,data1_size))));
        /* free memory */
        free(data1_);
        free(data2_);

## GSL function: double gsl_stats_covariance(const double data1[], const size_t stride1,const double data2[], const size_t stride2, const size_t n1);
void
stats_covariance(data1,data2)
        SV * data1;
        SV * data2;
    INIT:
        int i;
        AV * av;
        double * data1_;
        size_t data1_size;
        double * data2_;
        size_t data2_size;
    PPCODE:
        /* check valid ARRAYREF */
        data1_size = arr_ref_size(data1, "Warning: stats_covariance() - invalid 'data1' argument");
        data2_size = arr_ref_size(data2, "Warning: stats_covariance() - invalid 'data2' argument");
        if (data1_size<0 || data2_size<0) XSRETURN_UNDEF;
        /* check size */
        if (data1_size!=data2_size) {
          warn("Warning: stats_covariance() - array sizes differ: data1_size(%d), data2_size(%d)", data1_size, data2_size);
          XSRETURN_UNDEF;
        }
        /* allocate memory */
        data1_ = malloc(data1_size * sizeof(double));
        data2_ = malloc(data2_size * sizeof(double));
        if (!data1_ || !data2_) {
          warn("Warning: stats_covariance() - malloc failed");
          if (data1_) free(data1_);
          if (data2_) free(data2_);
          XSRETURN_UNDEF;
        }
        /* copy data */
        for(i=0, av=(AV*)SvRV(data1); i<data1_size; i++) data1_[i] = SvNV(*av_fetch(av,i,0));
        for(i=0, av=(AV*)SvRV(data2); i<data2_size; i++) data2_[i] = SvNV(*av_fetch(av,i,0));
        /* do the job */
        XPUSHs(sv_2mortal(newSVnv(gsl_stats_covariance(data1_,1,data2_,1,data1_size))));
        /* free memory */
        free(data1_);
        free(data2_);

## GSL function: double gsl_stats_covariance_m(const double data1[], const size_t stride1,const double data2[], const size_t stride2, const size_t n1, const double mean1, const double mean2);
void
stats_covariance_m(data1,data2,mean1,mean2)
        SV * data1;
        SV * data2;
        double mean1;
        double mean2;
    INIT:
        int i;
        AV * av;
        double * data1_;
        size_t data1_size;
        double * data2_;
        size_t data2_size;
    PPCODE:
        /* check valid ARRAYREF */
        data1_size = arr_ref_size(data1, "Warning: stats_covariance_m() - invalid 'data1' argument");
        data2_size = arr_ref_size(data2, "Warning: stats_covariance_m() - invalid 'data2' argument");
        if (data1_size<0 || data2_size<0) XSRETURN_UNDEF;
        /* check size */
        if (data1_size!=data2_size) {
          warn("Warning: stats_covariance_m() - array sizes differ: data1_size(%d), data2_size(%d)", data1_size, data2_size);
          XSRETURN_UNDEF;
        }
        /* allocate memory */
        data1_ = malloc(data1_size * sizeof(double));
        data2_ = malloc(data2_size * sizeof(double));
        if (!data1_ || !data2_) {
          warn("Warning: stats_covariance_m() - malloc failed");
          if (data1_) free(data1_);
          if (data2_) free(data2_);
          XSRETURN_UNDEF;
        }
        /* copy data */
        for(i=0, av=(AV*)SvRV(data1); i<data1_size; i++) data1_[i] = SvNV(*av_fetch(av,i,0));
        for(i=0, av=(AV*)SvRV(data2); i<data2_size; i++) data2_[i] = SvNV(*av_fetch(av,i,0));
        /* do the job */
        XPUSHs(sv_2mortal(newSVnv(gsl_stats_covariance_m(data1_,1,data2_,1,data1_size,mean1,mean2))));
        /* free memory */
        free(data1_);
        free(data2_);

## GSL function: double gsl_stats_kurtosis(const double data[], const size_t stride, const size_t n);
void
stats_kurtosis(data)
        SV * data;
    INIT:
        int i;
        AV * av;
        double * data_;
        size_t data_size;
    PPCODE:
        /* check valid ARRAYREF */
        data_size = arr_ref_size(data, "Warning: stats_kurtosis() - invalid 'data' argument");
        if (data_size<0) XSRETURN_UNDEF;
        /* allocate memory */
        data_ = malloc(data_size * sizeof(double));
        if (!data_) {
          warn("Warning: stats_kurtosis() - malloc failed");
          XSRETURN_UNDEF;
        }
        /* copy data */
        for(i=0, av=(AV*)SvRV(data); i<data_size; i++) data_[i] = SvNV(*av_fetch(av,i,0));
        /* do the job */
        XPUSHs(sv_2mortal(newSVnv(gsl_stats_kurtosis(data_,1,data_size))));
        /* free memory */
        free(data_);

## GSL function: double gsl_stats_kurtosis_m_sd(const double data[], const size_t stride, const size_t n, const double mean, const double sd);
void
stats_kurtosis_m_sd(data,mean,sd)
        SV * data;
        double mean;
        double sd;
    INIT:
        int i;
        AV * av;
        double * data_;
        size_t data_size;
    PPCODE:
        /* check valid ARRAYREF */
        data_size = arr_ref_size(data, "Warning: stats_kurtosis_m_sd() - invalid 'data' argument");
        if (data_size<0) XSRETURN_UNDEF;
        /* allocate memory */
        data_ = malloc(data_size * sizeof(double));
        if (!data_) {
          warn("Warning: stats_kurtosis_m_sd() - malloc failed");
          XSRETURN_UNDEF;
        }
        /* copy data */
        for(i=0, av=(AV*)SvRV(data); i<data_size; i++) data_[i] = SvNV(*av_fetch(av,i,0));
        /* do the job */
        XPUSHs(sv_2mortal(newSVnv(gsl_stats_kurtosis_m_sd(data_,1,data_size,mean,sd))));
        /* free memory */
        free(data_);

## GSL function: double gsl_stats_lag1_autocorrelation(const double data[], const size_t stride, const size_t n);
void
stats_lag1_autocorrelation(data)
        SV * data;
    INIT:
        int i;
        AV * av;
        double * data_;
        size_t data_size;
    PPCODE:
        /* check valid ARRAYREF */
        data_size = arr_ref_size(data, "Warning: stats_lag1_autocorrelation() - invalid 'data' argument");
        if (data_size<0) XSRETURN_UNDEF;
        /* allocate memory */
        data_ = malloc(data_size * sizeof(double));
        if (!data_) {
          warn("Warning: stats_lag1_autocorrelation() - malloc failed");
          XSRETURN_UNDEF;
        }
        /* copy data */
        for(i=0, av=(AV*)SvRV(data); i<data_size; i++) data_[i] = SvNV(*av_fetch(av,i,0));
        /* do the job */
        XPUSHs(sv_2mortal(newSVnv(gsl_stats_lag1_autocorrelation(data_,1,data_size))));
        /* free memory */
        free(data_);

## GSL function: double gsl_stats_lag1_autocorrelation_m(const double data[], const size_t stride, const size_t n, const double mean);
void
stats_lag1_autocorrelation_m(data,mean)
        SV * data;
        double mean;
    INIT:
        int i;
        AV * av;
        double * data_;
        size_t data_size;
    PPCODE:
        /* check valid ARRAYREF */
        data_size = arr_ref_size(data, "Warning: stats_lag1_autocorrelation_m() - invalid 'data' argument");
        if (data_size<0) XSRETURN_UNDEF;
        /* allocate memory */
        data_ = malloc(data_size * sizeof(double));
        if (!data_) {
          warn("Warning: stats_lag1_autocorrelation_m() - malloc failed");
          XSRETURN_UNDEF;
        }
        /* copy data */
        for(i=0, av=(AV*)SvRV(data); i<data_size; i++) data_[i] = SvNV(*av_fetch(av,i,0));
        /* do the job */
        XPUSHs(sv_2mortal(newSVnv(gsl_stats_lag1_autocorrelation_m(data_,1,data_size,mean))));
        /* free memory */
        free(data_);

## GSL function: double gsl_stats_max(const double data[], const size_t stride, const size_t n);
void
stats_max(data)
        SV * data;
    INIT:
        int i;
        AV * av;
        double * data_;
        size_t data_size;
    PPCODE:
        /* check valid ARRAYREF */
        data_size = arr_ref_size(data, "Warning: stats_max() - invalid 'data' argument");
        if (data_size<0) XSRETURN_UNDEF;
        /* allocate memory */
        data_ = malloc(data_size * sizeof(double));
        if (!data_) {
          warn("Warning: stats_max() - malloc failed");
          XSRETURN_UNDEF;
        }
        /* copy data */
        for(i=0, av=(AV*)SvRV(data); i<data_size; i++) data_[i] = SvNV(*av_fetch(av,i,0));
        /* do the job */
        XPUSHs(sv_2mortal(newSVnv(gsl_stats_max(data_,1,data_size))));
        /* free memory */
        free(data_);

## GSL function: size_t gsl_stats_max_index(const double data[], const size_t stride, const size_t n);
void
stats_max_index(data)
        SV * data;
    INIT:
        int i;
        AV * av;
        double * data_;
        size_t data_size;
    PPCODE:
        /* check valid ARRAYREF */
        data_size = arr_ref_size(data, "Warning: stats_max_index() - invalid 'data' argument");
        if (data_size<0) XSRETURN_UNDEF;
        /* allocate memory */
        data_ = malloc(data_size * sizeof(double));
        if (!data_) {
          warn("Warning: stats_max_index() - malloc failed");
          XSRETURN_UNDEF;
        }
        /* copy data */
        for(i=0, av=(AV*)SvRV(data); i<data_size; i++) data_[i] = SvNV(*av_fetch(av,i,0));
        /* do the job */
        XPUSHs(sv_2mortal(newSViv(gsl_stats_max_index(data_,1,data_size))));
        /* free memory */
        free(data_);

## GSL function: double gsl_stats_mean(const double data[], const size_t stride, const size_t n);
void
stats_mean(data)
        SV * data;
    INIT:
        int i;
        AV * av;
        double * data_;
        size_t data_size;
    PPCODE:
        /* check valid ARRAYREF */
        data_size = arr_ref_size(data, "Warning: stats_mean() - invalid 'data' argument");
        if (data_size<0) XSRETURN_UNDEF;
        /* allocate memory */
        data_ = malloc(data_size * sizeof(double));
        if (!data_) {
          warn("Warning: stats_mean() - malloc failed");
          XSRETURN_UNDEF;
        }
        /* copy data */
        for(i=0, av=(AV*)SvRV(data); i<data_size; i++) data_[i] = SvNV(*av_fetch(av,i,0));
        /* do the job */
        XPUSHs(sv_2mortal(newSVnv(gsl_stats_mean(data_,1,data_size))));
        /* free memory */
        free(data_);

## GSL function: double gsl_stats_min(const double data[], const size_t stride, const size_t n);
void
stats_min(data)
        SV * data;
    INIT:
        int i;
        AV * av;
        double * data_;
        size_t data_size;
    PPCODE:
        /* check valid ARRAYREF */
        data_size = arr_ref_size(data, "Warning: stats_min() - invalid 'data' argument");
        if (data_size<0) XSRETURN_UNDEF;
        /* allocate memory */
        data_ = malloc(data_size * sizeof(double));
        if (!data_) {
          warn("Warning: stats_min() - malloc failed");
          XSRETURN_UNDEF;
        }
        /* copy data */
        for(i=0, av=(AV*)SvRV(data); i<data_size; i++) data_[i] = SvNV(*av_fetch(av,i,0));
        /* do the job */
        XPUSHs(sv_2mortal(newSVnv(gsl_stats_min(data_,1,data_size))));
        /* free memory */
        free(data_);

## GSL function: size_t gsl_stats_min_index(const double data[], const size_t stride, const size_t n);
void
stats_min_index(data)
        SV * data;
    INIT:
        int i;
        AV * av;
        double * data_;
        size_t data_size;
    PPCODE:
        /* check valid ARRAYREF */
        data_size = arr_ref_size(data, "Warning: stats_min_index() - invalid 'data' argument");
        if (data_size<0) XSRETURN_UNDEF;
        /* allocate memory */
        data_ = malloc(data_size * sizeof(double));
        if (!data_) {
          warn("Warning: stats_min_index() - malloc failed");
          XSRETURN_UNDEF;
        }
        /* copy data */
        for(i=0, av=(AV*)SvRV(data); i<data_size; i++) data_[i] = SvNV(*av_fetch(av,i,0));
        /* do the job */
        XPUSHs(sv_2mortal(newSViv(gsl_stats_min_index(data_,1,data_size))));
        /* free memory */
        free(data_);

## GSL function: double gsl_stats_pvariance(const double data1[], const size_t stride1, const size_t n1, const double data2[], const size_t stride2, const size_t n2);
void
stats_pvariance(data1,data2)
        SV * data1;
        SV * data2;
    INIT:
        int i;
        AV * av;
        double * data1_;
        size_t data1_size;
        double * data2_;
        size_t data2_size;
    PPCODE:
        /* check valid ARRAYREF */
        data1_size = arr_ref_size(data1, "Warning: stats_pvariance() - invalid 'data1' argument");
        data2_size = arr_ref_size(data2, "Warning: stats_pvariance() - invalid 'data2' argument");
        if (data1_size<0 || data2_size<0) XSRETURN_UNDEF;
        /* check size */
        if (data1_size!=data2_size) {
          warn("Warning: stats_pvariance() - array sizes differ: data1_size(%d), data2_size(%d)", data1_size, data2_size);
          XSRETURN_UNDEF;
        }
        /* allocate memory */
        data1_ = malloc(data1_size * sizeof(double));
        data2_ = malloc(data2_size * sizeof(double));
        if (!data1_ || !data2_) {
          warn("Warning: stats_pvariance() - malloc failed");
          if (data1_) free(data1_);
          if (data2_) free(data2_);
          XSRETURN_UNDEF;
        }
        /* copy data */
        for(i=0, av=(AV*)SvRV(data1); i<data1_size; i++) data1_[i] = SvNV(*av_fetch(av,i,0));
        for(i=0, av=(AV*)SvRV(data2); i<data2_size; i++) data2_[i] = SvNV(*av_fetch(av,i,0));
        /* do the job */
        XPUSHs(sv_2mortal(newSVnv(gsl_stats_pvariance(data1_,1,data1_size,data2_,1,data2_size))));
        /* free memory */
        free(data1_);
        free(data2_);

## GSL function: double gsl_stats_sd(const double data[], const size_t stride, const size_t n);
void
stats_sd(data)
        SV * data;
    INIT:
        int i;
        AV * av;
        double * data_;
        size_t data_size;
    PPCODE:
        /* check valid ARRAYREF */
        data_size = arr_ref_size(data, "Warning: stats_sd() - invalid 'data' argument");
        if (data_size<0) XSRETURN_UNDEF;
        /* allocate memory */
        data_ = malloc(data_size * sizeof(double));
        if (!data_) {
          warn("Warning: stats_sd() - malloc failed");
          XSRETURN_UNDEF;
        }
        /* copy data */
        for(i=0, av=(AV*)SvRV(data); i<data_size; i++) data_[i] = SvNV(*av_fetch(av,i,0));
        /* do the job */
        XPUSHs(sv_2mortal(newSVnv(gsl_stats_sd(data_,1,data_size))));
        /* free memory */
        free(data_);

## GSL function: double gsl_stats_sd_m(const double data[], const size_t stride, const size_t n, const double mean);
void
stats_sd_m(data,mean)
        SV * data;
        double mean;
    INIT:
        int i;
        AV * av;
        double * data_;
        size_t data_size;
    PPCODE:
        /* check valid ARRAYREF */
        data_size = arr_ref_size(data, "Warning: stats_sd_m() - invalid 'data' argument");
        if (data_size<0) XSRETURN_UNDEF;
        /* allocate memory */
        data_ = malloc(data_size * sizeof(double));
        if (!data_) {
          warn("Warning: stats_sd_m() - malloc failed");
          XSRETURN_UNDEF;
        }
        /* copy data */
        for(i=0, av=(AV*)SvRV(data); i<data_size; i++) data_[i] = SvNV(*av_fetch(av,i,0));
        /* do the job */
        XPUSHs(sv_2mortal(newSVnv(gsl_stats_sd_m(data_,1,data_size,mean))));
        /* free memory */
        free(data_);

## GSL function: double gsl_stats_sd_with_fixed_mean(const double data[], const size_t stride, const size_t n, const double mean);
void
stats_sd_with_fixed_mean(data,mean)
        SV * data;
        double mean;
    INIT:
        int i;
        AV * av;
        double * data_;
        size_t data_size;
    PPCODE:
        /* check valid ARRAYREF */
        data_size = arr_ref_size(data, "Warning: stats_sd_with_fixed_mean() - invalid 'data' argument");
        if (data_size<0) XSRETURN_UNDEF;
        /* allocate memory */
        data_ = malloc(data_size * sizeof(double));
        if (!data_) {
          warn("Warning: stats_sd_with_fixed_mean() - malloc failed");
          XSRETURN_UNDEF;
        }
        /* copy data */
        for(i=0, av=(AV*)SvRV(data); i<data_size; i++) data_[i] = SvNV(*av_fetch(av,i,0));
        /* do the job */
        XPUSHs(sv_2mortal(newSVnv(gsl_stats_sd_with_fixed_mean(data_,1,data_size,mean))));
        /* free memory */
        free(data_);

## GSL function: double gsl_stats_skew(const double data[], const size_t stride, const size_t n);
void
stats_skew(data)
        SV * data;
    INIT:
        int i;
        AV * av;
        double * data_;
        size_t data_size;
    PPCODE:
        /* check valid ARRAYREF */
        data_size = arr_ref_size(data, "Warning: stats_skew() - invalid 'data' argument");
        if (data_size<0) XSRETURN_UNDEF;
        /* allocate memory */
        data_ = malloc(data_size * sizeof(double));
        if (!data_) {
          warn("Warning: stats_skew() - malloc failed");
          XSRETURN_UNDEF;
        }
        /* copy data */
        for(i=0, av=(AV*)SvRV(data); i<data_size; i++) data_[i] = SvNV(*av_fetch(av,i,0));
        /* do the job */
        XPUSHs(sv_2mortal(newSVnv(gsl_stats_skew(data_,1,data_size))));
        /* free memory */
        free(data_);

## GSL function: double gsl_stats_skew_m_sd(const double data[], const size_t stride, const size_t n, const double mean, const double sd);
void
stats_skew_m_sd(data,mean,sd)
        SV * data;
        double mean;
        double sd;
    INIT:
        int i;
        AV * av;
        double * data_;
        size_t data_size;
    PPCODE:
        /* check valid ARRAYREF */
        data_size = arr_ref_size(data, "Warning: stats_skew_m_sd() - invalid 'data' argument");
        if (data_size<0) XSRETURN_UNDEF;
        /* allocate memory */
        data_ = malloc(data_size * sizeof(double));
        if (!data_) {
          warn("Warning: stats_skew_m_sd() - malloc failed");
          XSRETURN_UNDEF;
        }
        /* copy data */
        for(i=0, av=(AV*)SvRV(data); i<data_size; i++) data_[i] = SvNV(*av_fetch(av,i,0));
        /* do the job */
        XPUSHs(sv_2mortal(newSVnv(gsl_stats_skew_m_sd(data_,1,data_size,mean,sd))));
        /* free memory */
        free(data_);

## GSL function: double gsl_stats_tss(const double data[], const size_t stride, const size_t n);
void
stats_tss(data)
        SV * data;
    INIT:
        int i;
        AV * av;
        double * data_;
        size_t data_size;
    PPCODE:
        /* check valid ARRAYREF */
        data_size = arr_ref_size(data, "Warning: stats_tss() - invalid 'data' argument");
        if (data_size<0) XSRETURN_UNDEF;
        /* allocate memory */
        data_ = malloc(data_size * sizeof(double));
        if (!data_) {
          warn("Warning: stats_tss() - malloc failed");
          XSRETURN_UNDEF;
        }
        /* copy data */
        for(i=0, av=(AV*)SvRV(data); i<data_size; i++) data_[i] = SvNV(*av_fetch(av,i,0));
        /* do the job */
        XPUSHs(sv_2mortal(newSVnv(gsl_stats_tss(data_,1,data_size))));
        /* free memory */
        free(data_);

## GSL function: double gsl_stats_tss_m(const double data[], const size_t stride, const size_t n, const double mean);
void
stats_tss_m(data,mean)
        SV * data;
        double mean;
    INIT:
        int i;
        AV * av;
        double * data_;
        size_t data_size;
    PPCODE:
        /* check valid ARRAYREF */
        data_size = arr_ref_size(data, "Warning: stats_tss_m() - invalid 'data' argument");
        if (data_size<0) XSRETURN_UNDEF;
        /* allocate memory */
        data_ = malloc(data_size * sizeof(double));
        if (!data_) {
          warn("Warning: stats_tss_m() - malloc failed");
          XSRETURN_UNDEF;
        }
        /* copy data */
        for(i=0, av=(AV*)SvRV(data); i<data_size; i++) data_[i] = SvNV(*av_fetch(av,i,0));
        /* do the job */
        XPUSHs(sv_2mortal(newSVnv(gsl_stats_tss_m(data_,1,data_size,mean))));
        /* free memory */
        free(data_);

## GSL function: double gsl_stats_ttest(const double data1[], const size_t stride1, const size_t n1, const double data2[], const size_t stride2, const size_t n2);
void
stats_ttest(data1,data2)
        SV * data1;
        SV * data2;
    INIT:
        int i;
        AV * av;
        double * data1_;
        size_t data1_size;
        double * data2_;
        size_t data2_size;
    PPCODE:
        /* check valid ARRAYREF */
        data1_size = arr_ref_size(data1, "Warning: stats_ttest() - invalid 'data1' argument");
        data2_size = arr_ref_size(data2, "Warning: stats_ttest() - invalid 'data2' argument");
        if (data1_size<0 || data2_size<0) XSRETURN_UNDEF;
        /* check size */
        if (data1_size!=data2_size) {
          warn("Warning: stats_ttest() - array sizes differ: data1_size(%d), data2_size(%d)", data1_size, data2_size);
          XSRETURN_UNDEF;
        }
        /* allocate memory */
        data1_ = malloc(data1_size * sizeof(double));
        data2_ = malloc(data2_size * sizeof(double));
        if (!data1_ || !data2_) {
          warn("Warning: stats_ttest() - malloc failed");
          if (data1_) free(data1_);
          if (data2_) free(data2_);
          XSRETURN_UNDEF;
        }
        /* copy data */
        for(i=0, av=(AV*)SvRV(data1); i<data1_size; i++) data1_[i] = SvNV(*av_fetch(av,i,0));
        for(i=0, av=(AV*)SvRV(data2); i<data2_size; i++) data2_[i] = SvNV(*av_fetch(av,i,0));
        /* do the job */
        XPUSHs(sv_2mortal(newSVnv(gsl_stats_ttest(data1_,1,data1_size,data2_,1,data2_size))));
        /* free memory */
        free(data1_);
        free(data2_);

## GSL function: double gsl_stats_variance(const double data[], const size_t stride, const size_t n);
void
stats_variance(data)
        SV * data;
    INIT:
        int i;
        AV * av;
        double * data_;
        size_t data_size;
    PPCODE:
        /* check valid ARRAYREF */
        data_size = arr_ref_size(data, "Warning: stats_variance() - invalid 'data' argument");
        if (data_size<0) XSRETURN_UNDEF;
        /* allocate memory */
        data_ = malloc(data_size * sizeof(double));
        if (!data_) {
          warn("Warning: stats_variance() - malloc failed");
          XSRETURN_UNDEF;
        }
        /* copy data */
        for(i=0, av=(AV*)SvRV(data); i<data_size; i++) data_[i] = SvNV(*av_fetch(av,i,0));
        /* do the job */
        XPUSHs(sv_2mortal(newSVnv(gsl_stats_variance(data_,1,data_size))));
        /* free memory */
        free(data_);

## GSL function: double gsl_stats_variance_m(const double data[], const size_t stride, const size_t n, const double mean);
void
stats_variance_m(data,mean)
        SV * data;
        double mean;
    INIT:
        int i;
        AV * av;
        double * data_;
        size_t data_size;
    PPCODE:
        /* check valid ARRAYREF */
        data_size = arr_ref_size(data, "Warning: stats_variance_m() - invalid 'data' argument");
        if (data_size<0) XSRETURN_UNDEF;
        /* allocate memory */
        data_ = malloc(data_size * sizeof(double));
        if (!data_) {
          warn("Warning: stats_variance_m() - malloc failed");
          XSRETURN_UNDEF;
        }
        /* copy data */
        for(i=0, av=(AV*)SvRV(data); i<data_size; i++) data_[i] = SvNV(*av_fetch(av,i,0));
        /* do the job */
        XPUSHs(sv_2mortal(newSVnv(gsl_stats_variance_m(data_,1,data_size,mean))));
        /* free memory */
        free(data_);

## GSL function: double gsl_stats_variance_with_fixed_mean(const double data[], const size_t stride, const size_t n, const double mean);
void
stats_variance_with_fixed_mean(data,mean)
        SV * data;
        double mean;
    INIT:
        int i;
        AV * av;
        double * data_;
        size_t data_size;
    PPCODE:
        /* check valid ARRAYREF */
        data_size = arr_ref_size(data, "Warning: stats_variance_with_fixed_mean() - invalid 'data' argument");
        if (data_size<0) XSRETURN_UNDEF;
        /* allocate memory */
        data_ = malloc(data_size * sizeof(double));
        if (!data_) {
          warn("Warning: stats_variance_with_fixed_mean() - malloc failed");
          XSRETURN_UNDEF;
        }
        /* copy data */
        for(i=0, av=(AV*)SvRV(data); i<data_size; i++) data_[i] = SvNV(*av_fetch(av,i,0));
        /* do the job */
        XPUSHs(sv_2mortal(newSVnv(gsl_stats_variance_with_fixed_mean(data_,1,data_size,mean))));
        /* free memory */
        free(data_);

## GSL function: double gsl_stats_wabsdev(const double w[], const size_t wstride, const double data[], const size_t stride, const size_t n);
void
stats_wabsdev(w,data)
        SV * w;
        SV * data;
    INIT:
        int i;
        AV * av;
        double * w_;
        size_t w_size;
        double * data_;
        size_t data_size;
    PPCODE:
        /* check valid ARRAYREF */
        w_size = arr_ref_size(w, "Warning: stats_wabsdev() - invalid 'w' argument");
        data_size = arr_ref_size(data, "Warning: stats_wabsdev() - invalid 'data' argument");
        if (w_size<0 || data_size<0) XSRETURN_UNDEF;
        /* check size */
        if (w_size!=data_size) {
          warn("Warning: stats_wabsdev() - array sizes differ: w_size(%d), data_size(%d)", w_size, data_size);
          XSRETURN_UNDEF;
        }
        /* allocate memory */
        w_ = malloc(w_size * sizeof(double));
        data_ = malloc(data_size * sizeof(double));
        if (!w_ || !data_) {
          warn("Warning: stats_wabsdev() - malloc failed");
          if (w_) free(w_);
          if (data_) free(data_);
          XSRETURN_UNDEF;
        }
        /* copy data */
        for(i=0, av=(AV*)SvRV(w); i<w_size; i++) w_[i] = SvNV(*av_fetch(av,i,0));
        for(i=0, av=(AV*)SvRV(data); i<data_size; i++) data_[i] = SvNV(*av_fetch(av,i,0));
        /* do the job */
        XPUSHs(sv_2mortal(newSVnv(gsl_stats_wabsdev(w_,1,data_,1,data_size))));
        /* free memory */
        free(w_);
        free(data_);

## GSL function: double gsl_stats_wabsdev_m(const double w[], const size_t wstride, const double data[], const size_t stride, const size_t n, const double wmean);
void
stats_wabsdev_m(w,data,wmean)
        SV * w;
        SV * data;
        double wmean;
    INIT:
        int i;
        AV * av;
        double * w_;
        size_t w_size;
        double * data_;
        size_t data_size;
    PPCODE:
        /* check valid ARRAYREF */
        w_size = arr_ref_size(w, "Warning: stats_wabsdev_m() - invalid 'w' argument");
        data_size = arr_ref_size(data, "Warning: stats_wabsdev_m() - invalid 'data' argument");
        if (w_size<0 || data_size<0) XSRETURN_UNDEF;
        /* check size */
        if (w_size!=data_size) {
          warn("Warning: stats_wabsdev_m() - array sizes differ: w_size(%d), data_size(%d)", w_size, data_size);
          XSRETURN_UNDEF;
        }
        /* allocate memory */
        w_ = malloc(w_size * sizeof(double));
        data_ = malloc(data_size * sizeof(double));
        if (!w_ || !data_) {
          warn("Warning: stats_wabsdev_m() - malloc failed");
          if (w_) free(w_);
          if (data_) free(data_);
          XSRETURN_UNDEF;
        }
        /* copy data */
        for(i=0, av=(AV*)SvRV(w); i<w_size; i++) w_[i] = SvNV(*av_fetch(av,i,0));
        for(i=0, av=(AV*)SvRV(data); i<data_size; i++) data_[i] = SvNV(*av_fetch(av,i,0));
        /* do the job */
        XPUSHs(sv_2mortal(newSVnv(gsl_stats_wabsdev_m(w_,1,data_,1,data_size,wmean))));
        /* free memory */
        free(w_);
        free(data_);

## GSL function: double gsl_stats_wkurtosis(const double w[], const size_t wstride, const double data[], const size_t stride, const size_t n);
void
stats_wkurtosis(w,data)
        SV * w;
        SV * data;
    INIT:
        int i;
        AV * av;
        double * w_;
        size_t w_size;
        double * data_;
        size_t data_size;
    PPCODE:
        /* check valid ARRAYREF */
        w_size = arr_ref_size(w, "Warning: stats_wkurtosis() - invalid 'w' argument");
        data_size = arr_ref_size(data, "Warning: stats_wkurtosis() - invalid 'data' argument");
        if (w_size<0 || data_size<0) XSRETURN_UNDEF;
        /* check size */
        if (w_size!=data_size) {
          warn("Warning: stats_wkurtosis() - array sizes differ: w_size(%d), data_size(%d)", w_size, data_size);
          XSRETURN_UNDEF;
        }
        /* allocate memory */
        w_ = malloc(w_size * sizeof(double));
        data_ = malloc(data_size * sizeof(double));
        if (!w_ || !data_) {
          warn("Warning: stats_wkurtosis() - malloc failed");
          if (w_) free(w_);
          if (data_) free(data_);
          XSRETURN_UNDEF;
        }
        /* copy data */
        for(i=0, av=(AV*)SvRV(w); i<w_size; i++) w_[i] = SvNV(*av_fetch(av,i,0));
        for(i=0, av=(AV*)SvRV(data); i<data_size; i++) data_[i] = SvNV(*av_fetch(av,i,0));
        /* do the job */
        XPUSHs(sv_2mortal(newSVnv(gsl_stats_wkurtosis(w_,1,data_,1,data_size))));
        /* free memory */
        free(w_);
        free(data_);

## GSL function: double gsl_stats_wkurtosis_m_sd(const double w[], const size_t wstride, const double data[], const size_t stride, const size_t n, const double wmean, const double wsd);
void
stats_wkurtosis_m_sd(w,data,wmean,wsd)
        SV * w;
        SV * data;
        double wmean;
        double wsd;
    INIT:
        int i;
        AV * av;
        double * w_;
        size_t w_size;
        double * data_;
        size_t data_size;
    PPCODE:
        /* check valid ARRAYREF */
        w_size = arr_ref_size(w, "Warning: stats_wkurtosis_m_sd() - invalid 'w' argument");
        data_size = arr_ref_size(data, "Warning: stats_wkurtosis_m_sd() - invalid 'data' argument");
        if (w_size<0 || data_size<0) XSRETURN_UNDEF;
        /* check size */
        if (w_size!=data_size) {
          warn("Warning: stats_wkurtosis_m_sd() - array sizes differ: w_size(%d), data_size(%d)", w_size, data_size);
          XSRETURN_UNDEF;
        }
        /* allocate memory */
        w_ = malloc(w_size * sizeof(double));
        data_ = malloc(data_size * sizeof(double));
        if (!w_ || !data_) {
          warn("Warning: stats_wkurtosis_m_sd() - malloc failed");
          if (w_) free(w_);
          if (data_) free(data_);
          XSRETURN_UNDEF;
        }
        /* copy data */
        for(i=0, av=(AV*)SvRV(w); i<w_size; i++) w_[i] = SvNV(*av_fetch(av,i,0));
        for(i=0, av=(AV*)SvRV(data); i<data_size; i++) data_[i] = SvNV(*av_fetch(av,i,0));
        /* do the job */
        XPUSHs(sv_2mortal(newSVnv(gsl_stats_wkurtosis_m_sd(w_,1,data_,1,data_size,wmean,wsd))));
        /* free memory */
        free(w_);
        free(data_);

## GSL function: double gsl_stats_wmean(const double w[], const size_t wstride, const double data[], const size_t stride, const size_t n);
void
stats_wmean(w,data)
        SV * w;
        SV * data;
    INIT:
        int i;
        AV * av;
        double * w_;
        size_t w_size;
        double * data_;
        size_t data_size;
    PPCODE:
        /* check valid ARRAYREF */
        w_size = arr_ref_size(w, "Warning: stats_wmean() - invalid 'w' argument");
        data_size = arr_ref_size(data, "Warning: stats_wmean() - invalid 'data' argument");
        if (w_size<0 || data_size<0) XSRETURN_UNDEF;
        /* check size */
        if (w_size!=data_size) {
          warn("Warning: stats_wmean() - array sizes differ: w_size(%d), data_size(%d)", w_size, data_size);
          XSRETURN_UNDEF;
        }
        /* allocate memory */
        w_ = malloc(w_size * sizeof(double));
        data_ = malloc(data_size * sizeof(double));
        if (!w_ || !data_) {
          warn("Warning: stats_wmean() - malloc failed");
          if (w_) free(w_);
          if (data_) free(data_);
          XSRETURN_UNDEF;
        }
        /* copy data */
        for(i=0, av=(AV*)SvRV(w); i<w_size; i++) w_[i] = SvNV(*av_fetch(av,i,0));
        for(i=0, av=(AV*)SvRV(data); i<data_size; i++) data_[i] = SvNV(*av_fetch(av,i,0));
        /* do the job */
        XPUSHs(sv_2mortal(newSVnv(gsl_stats_wmean(w_,1,data_,1,data_size))));
        /* free memory */
        free(w_);
        free(data_);

## GSL function: double gsl_stats_wsd(const double w[], const size_t wstride, const double data[], const size_t stride, const size_t n);
void
stats_wsd(w,data)
        SV * w;
        SV * data;
    INIT:
        int i;
        AV * av;
        double * w_;
        size_t w_size;
        double * data_;
        size_t data_size;
    PPCODE:
        /* check valid ARRAYREF */
        w_size = arr_ref_size(w, "Warning: stats_wsd() - invalid 'w' argument");
        data_size = arr_ref_size(data, "Warning: stats_wsd() - invalid 'data' argument");
        if (w_size<0 || data_size<0) XSRETURN_UNDEF;
        /* check size */
        if (w_size!=data_size) {
          warn("Warning: stats_wsd() - array sizes differ: w_size(%d), data_size(%d)", w_size, data_size);
          XSRETURN_UNDEF;
        }
        /* allocate memory */
        w_ = malloc(w_size * sizeof(double));
        data_ = malloc(data_size * sizeof(double));
        if (!w_ || !data_) {
          warn("Warning: stats_wsd() - malloc failed");
          if (w_) free(w_);
          if (data_) free(data_);
          XSRETURN_UNDEF;
        }
        /* copy data */
        for(i=0, av=(AV*)SvRV(w); i<w_size; i++) w_[i] = SvNV(*av_fetch(av,i,0));
        for(i=0, av=(AV*)SvRV(data); i<data_size; i++) data_[i] = SvNV(*av_fetch(av,i,0));
        /* do the job */
        XPUSHs(sv_2mortal(newSVnv(gsl_stats_wsd(w_,1,data_,1,data_size))));
        /* free memory */
        free(w_);
        free(data_);

## GSL function: double gsl_stats_wsd_m(const double w[], const size_t wstride, const double data[], const size_t stride, const size_t n, const double wmean);
void
stats_wsd_m(w,data,wmean)
        SV * w;
        SV * data;
        double wmean;
    INIT:
        int i;
        AV * av;
        double * w_;
        size_t w_size;
        double * data_;
        size_t data_size;
    PPCODE:
        /* check valid ARRAYREF */
        w_size = arr_ref_size(w, "Warning: stats_wsd_m() - invalid 'w' argument");
        data_size = arr_ref_size(data, "Warning: stats_wsd_m() - invalid 'data' argument");
        if (w_size<0 || data_size<0) XSRETURN_UNDEF;
        /* check size */
        if (w_size!=data_size) {
          warn("Warning: stats_wsd_m() - array sizes differ: w_size(%d), data_size(%d)", w_size, data_size);
          XSRETURN_UNDEF;
        }
        /* allocate memory */
        w_ = malloc(w_size * sizeof(double));
        data_ = malloc(data_size * sizeof(double));
        if (!w_ || !data_) {
          warn("Warning: stats_wsd_m() - malloc failed");
          if (w_) free(w_);
          if (data_) free(data_);
          XSRETURN_UNDEF;
        }
        /* copy data */
        for(i=0, av=(AV*)SvRV(w); i<w_size; i++) w_[i] = SvNV(*av_fetch(av,i,0));
        for(i=0, av=(AV*)SvRV(data); i<data_size; i++) data_[i] = SvNV(*av_fetch(av,i,0));
        /* do the job */
        XPUSHs(sv_2mortal(newSVnv(gsl_stats_wsd_m(w_,1,data_,1,data_size,wmean))));
        /* free memory */
        free(w_);
        free(data_);

## GSL function: double gsl_stats_wsd_with_fixed_mean(const double w[], const size_t wstride, const double data[], const size_t stride, const size_t n, const double mean);
void
stats_wsd_with_fixed_mean(w,data,mean)
        SV * w;
        SV * data;
        double mean;
    INIT:
        int i;
        AV * av;
        double * w_;
        size_t w_size;
        double * data_;
        size_t data_size;
    PPCODE:
        /* check valid ARRAYREF */
        w_size = arr_ref_size(w, "Warning: stats_wsd_with_fixed_mean() - invalid 'w' argument");
        data_size = arr_ref_size(data, "Warning: stats_wsd_with_fixed_mean() - invalid 'data' argument");
        if (w_size<0 || data_size<0) XSRETURN_UNDEF;
        /* check size */
        if (w_size!=data_size) {
          warn("Warning: stats_wsd_with_fixed_mean() - array sizes differ: w_size(%d), data_size(%d)", w_size, data_size);
          XSRETURN_UNDEF;
        }
        /* allocate memory */
        w_ = malloc(w_size * sizeof(double));
        data_ = malloc(data_size * sizeof(double));
        if (!w_ || !data_) {
          warn("Warning: stats_wsd_with_fixed_mean() - malloc failed");
          if (w_) free(w_);
          if (data_) free(data_);
          XSRETURN_UNDEF;
        }
        /* copy data */
        for(i=0, av=(AV*)SvRV(w); i<w_size; i++) w_[i] = SvNV(*av_fetch(av,i,0));
        for(i=0, av=(AV*)SvRV(data); i<data_size; i++) data_[i] = SvNV(*av_fetch(av,i,0));
        /* do the job */
        XPUSHs(sv_2mortal(newSVnv(gsl_stats_wsd_with_fixed_mean(w_,1,data_,1,data_size,mean))));
        /* free memory */
        free(w_);
        free(data_);

## GSL function: double gsl_stats_wskew(const double w[], const size_t wstride, const double data[], const size_t stride, const size_t n);
void
stats_wskew(w,data)
        SV * w;
        SV * data;
    INIT:
        int i;
        AV * av;
        double * w_;
        size_t w_size;
        double * data_;
        size_t data_size;
    PPCODE:
        /* check valid ARRAYREF */
        w_size = arr_ref_size(w, "Warning: stats_wskew() - invalid 'w' argument");
        data_size = arr_ref_size(data, "Warning: stats_wskew() - invalid 'data' argument");
        if (w_size<0 || data_size<0) XSRETURN_UNDEF;
        /* check size */
        if (w_size!=data_size) {
          warn("Warning: stats_wskew() - array sizes differ: w_size(%d), data_size(%d)", w_size, data_size);
          XSRETURN_UNDEF;
        }
        /* allocate memory */
        w_ = malloc(w_size * sizeof(double));
        data_ = malloc(data_size * sizeof(double));
        if (!w_ || !data_) {
          warn("Warning: stats_wskew() - malloc failed");
          if (w_) free(w_);
          if (data_) free(data_);
          XSRETURN_UNDEF;
        }
        /* copy data */
        for(i=0, av=(AV*)SvRV(w); i<w_size; i++) w_[i] = SvNV(*av_fetch(av,i,0));
        for(i=0, av=(AV*)SvRV(data); i<data_size; i++) data_[i] = SvNV(*av_fetch(av,i,0));
        /* do the job */
        XPUSHs(sv_2mortal(newSVnv(gsl_stats_wskew(w_,1,data_,1,data_size))));
        /* free memory */
        free(w_);
        free(data_);

## GSL function: double gsl_stats_wskew_m_sd(const double w[], const size_t wstride, const double data[], const size_t stride, const size_t n, const double wmean, const double wsd);
void
stats_wskew_m_sd(w,data,wmean,wsd)
        SV * w;
        SV * data;
        double wmean;
        double wsd;
    INIT:
        int i;
        AV * av;
        double * w_;
        size_t w_size;
        double * data_;
        size_t data_size;
    PPCODE:
        /* check valid ARRAYREF */
        w_size = arr_ref_size(w, "Warning: stats_wskew_m_sd() - invalid 'w' argument");
        data_size = arr_ref_size(data, "Warning: stats_wskew_m_sd() - invalid 'data' argument");
        if (w_size<0 || data_size<0) XSRETURN_UNDEF;
        /* check size */
        if (w_size!=data_size) {
          warn("Warning: stats_wskew_m_sd() - array sizes differ: w_size(%d), data_size(%d)", w_size, data_size);
          XSRETURN_UNDEF;
        }
        /* allocate memory */
        w_ = malloc(w_size * sizeof(double));
        data_ = malloc(data_size * sizeof(double));
        if (!w_ || !data_) {
          warn("Warning: stats_wskew_m_sd() - malloc failed");
          if (w_) free(w_);
          if (data_) free(data_);
          XSRETURN_UNDEF;
        }
        /* copy data */
        for(i=0, av=(AV*)SvRV(w); i<w_size; i++) w_[i] = SvNV(*av_fetch(av,i,0));
        for(i=0, av=(AV*)SvRV(data); i<data_size; i++) data_[i] = SvNV(*av_fetch(av,i,0));
        /* do the job */
        XPUSHs(sv_2mortal(newSVnv(gsl_stats_wskew_m_sd(w_,1,data_,1,data_size,wmean,wsd))));
        /* free memory */
        free(w_);
        free(data_);

## GSL function: double gsl_stats_wtss(const double w[], const size_t wstride, const double data[], const size_t stride, const size_t n);
void
stats_wtss(w,data)
        SV * w;
        SV * data;
    INIT:
        int i;
        AV * av;
        double * w_;
        size_t w_size;
        double * data_;
        size_t data_size;
    PPCODE:
        /* check valid ARRAYREF */
        w_size = arr_ref_size(w, "Warning: stats_wtss() - invalid 'w' argument");
        data_size = arr_ref_size(data, "Warning: stats_wtss() - invalid 'data' argument");
        if (w_size<0 || data_size<0) XSRETURN_UNDEF;
        /* check size */
        if (w_size!=data_size) {
          warn("Warning: stats_wtss() - array sizes differ: w_size(%d), data_size(%d)", w_size, data_size);
          XSRETURN_UNDEF;
        }
        /* allocate memory */
        w_ = malloc(w_size * sizeof(double));
        data_ = malloc(data_size * sizeof(double));
        if (!w_ || !data_) {
          warn("Warning: stats_wtss() - malloc failed");
          if (w_) free(w_);
          if (data_) free(data_);
          XSRETURN_UNDEF;
        }
        /* copy data */
        for(i=0, av=(AV*)SvRV(w); i<w_size; i++) w_[i] = SvNV(*av_fetch(av,i,0));
        for(i=0, av=(AV*)SvRV(data); i<data_size; i++) data_[i] = SvNV(*av_fetch(av,i,0));
        /* do the job */
        XPUSHs(sv_2mortal(newSVnv(gsl_stats_wtss(w_,1,data_,1,data_size))));
        /* free memory */
        free(w_);
        free(data_);

## GSL function: double gsl_stats_wtss_m(const double w[], const size_t wstride, const double data[], const size_t stride, const size_t n, const double wmean);
void
stats_wtss_m(w,data,wmean)
        SV * w;
        SV * data;
        double wmean;
    INIT:
        int i;
        AV * av;
        double * w_;
        size_t w_size;
        double * data_;
        size_t data_size;
    PPCODE:
        /* check valid ARRAYREF */
        w_size = arr_ref_size(w, "Warning: stats_wtss_m() - invalid 'w' argument");
        data_size = arr_ref_size(data, "Warning: stats_wtss_m() - invalid 'data' argument");
        if (w_size<0 || data_size<0) XSRETURN_UNDEF;
        /* check size */
        if (w_size!=data_size) {
          warn("Warning: stats_wtss_m() - array sizes differ: w_size(%d), data_size(%d)", w_size, data_size);
          XSRETURN_UNDEF;
        }
        /* allocate memory */
        w_ = malloc(w_size * sizeof(double));
        data_ = malloc(data_size * sizeof(double));
        if (!w_ || !data_) {
          warn("Warning: stats_wtss_m() - malloc failed");
          if (w_) free(w_);
          if (data_) free(data_);
          XSRETURN_UNDEF;
        }
        /* copy data */
        for(i=0, av=(AV*)SvRV(w); i<w_size; i++) w_[i] = SvNV(*av_fetch(av,i,0));
        for(i=0, av=(AV*)SvRV(data); i<data_size; i++) data_[i] = SvNV(*av_fetch(av,i,0));
        /* do the job */
        XPUSHs(sv_2mortal(newSVnv(gsl_stats_wtss_m(w_,1,data_,1,data_size,wmean))));
        /* free memory */
        free(w_);
        free(data_);

## GSL function: double gsl_stats_wvariance(const double w[], const size_t wstride, const double data[], const size_t stride, const size_t n);
void
stats_wvariance(w,data)
        SV * w;
        SV * data;
    INIT:
        int i;
        AV * av;
        double * w_;
        size_t w_size;
        double * data_;
        size_t data_size;
    PPCODE:
        /* check valid ARRAYREF */
        w_size = arr_ref_size(w, "Warning: stats_wvariance() - invalid 'w' argument");
        data_size = arr_ref_size(data, "Warning: stats_wvariance() - invalid 'data' argument");
        if (w_size<0 || data_size<0) XSRETURN_UNDEF;
        /* check size */
        if (w_size!=data_size) {
          warn("Warning: stats_wvariance() - array sizes differ: w_size(%d), data_size(%d)", w_size, data_size);
          XSRETURN_UNDEF;
        }
        /* allocate memory */
        w_ = malloc(w_size * sizeof(double));
        data_ = malloc(data_size * sizeof(double));
        if (!w_ || !data_) {
          warn("Warning: stats_wvariance() - malloc failed");
          if (w_) free(w_);
          if (data_) free(data_);
          XSRETURN_UNDEF;
        }
        /* copy data */
        for(i=0, av=(AV*)SvRV(w); i<w_size; i++) w_[i] = SvNV(*av_fetch(av,i,0));
        for(i=0, av=(AV*)SvRV(data); i<data_size; i++) data_[i] = SvNV(*av_fetch(av,i,0));
        /* do the job */
        XPUSHs(sv_2mortal(newSVnv(gsl_stats_wvariance(w_,1,data_,1,data_size))));
        /* free memory */
        free(w_);
        free(data_);

## GSL function: double gsl_stats_wvariance_m(const double w[], const size_t wstride, const double data[], const size_t stride, const size_t n, const double wmean);
void
stats_wvariance_m(w,data,wmean)
        SV * w;
        SV * data;
        double wmean;
    INIT:
        int i;
        AV * av;
        double * w_;
        size_t w_size;
        double * data_;
        size_t data_size;
    PPCODE:
        /* check valid ARRAYREF */
        w_size = arr_ref_size(w, "Warning: stats_wvariance_m() - invalid 'w' argument");
        data_size = arr_ref_size(data, "Warning: stats_wvariance_m() - invalid 'data' argument");
        if (w_size<0 || data_size<0) XSRETURN_UNDEF;
        /* check size */
        if (w_size!=data_size) {
          warn("Warning: stats_wvariance_m() - array sizes differ: w_size(%d), data_size(%d)", w_size, data_size);
          XSRETURN_UNDEF;
        }
        /* allocate memory */
        w_ = malloc(w_size * sizeof(double));
        data_ = malloc(data_size * sizeof(double));
        if (!w_ || !data_) {
          warn("Warning: stats_wvariance_m() - malloc failed");
          if (w_) free(w_);
          if (data_) free(data_);
          XSRETURN_UNDEF;
        }
        /* copy data */
        for(i=0, av=(AV*)SvRV(w); i<w_size; i++) w_[i] = SvNV(*av_fetch(av,i,0));
        for(i=0, av=(AV*)SvRV(data); i<data_size; i++) data_[i] = SvNV(*av_fetch(av,i,0));
        /* do the job */
        XPUSHs(sv_2mortal(newSVnv(gsl_stats_wvariance_m(w_,1,data_,1,data_size,wmean))));
        /* free memory */
        free(w_);
        free(data_);

## GSL function: double gsl_stats_wvariance_with_fixed_mean(const double w[], const size_t wstride, const double data[], const size_t stride, const size_t n, const double mean);
void
stats_wvariance_with_fixed_mean(w,data,mean)
        SV * w;
        SV * data;
        double mean;
    INIT:
        int i;
        AV * av;
        double * w_;
        size_t w_size;
        double * data_;
        size_t data_size;
    PPCODE:
        /* check valid ARRAYREF */
        w_size = arr_ref_size(w, "Warning: stats_wvariance_with_fixed_mean() - invalid 'w' argument");
        data_size = arr_ref_size(data, "Warning: stats_wvariance_with_fixed_mean() - invalid 'data' argument");
        if (w_size<0 || data_size<0) XSRETURN_UNDEF;
        /* check size */
        if (w_size!=data_size) {
          warn("Warning: stats_wvariance_with_fixed_mean() - array sizes differ: w_size(%d), data_size(%d)", w_size, data_size);
          XSRETURN_UNDEF;
        }
        /* allocate memory */
        w_ = malloc(w_size * sizeof(double));
        data_ = malloc(data_size * sizeof(double));
        if (!w_ || !data_) {
          warn("Warning: stats_wvariance_with_fixed_mean() - malloc failed");
          if (w_) free(w_);
          if (data_) free(data_);
          XSRETURN_UNDEF;
        }
        /* copy data */
        for(i=0, av=(AV*)SvRV(w); i<w_size; i++) w_[i] = SvNV(*av_fetch(av,i,0));
        for(i=0, av=(AV*)SvRV(data); i<data_size; i++) data_[i] = SvNV(*av_fetch(av,i,0));
        /* do the job */
        XPUSHs(sv_2mortal(newSVnv(gsl_stats_wvariance_with_fixed_mean(w_,1,data_,1,data_size,mean))));
        /* free memory */
        free(w_);
        free(data_);

###### generated part - end ######

