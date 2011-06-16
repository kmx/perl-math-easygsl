package Math::EasyGSL::Statistics;

@ISA = qw/ DynaLoader /;
require DynaLoader;

bootstrap Math::EasyGSL::Statistics;

use Exporter 'import'; # gives you Exporter's import() method directly
@EXPORT_OK = qw(
stats_median
stats_minmax
stats_minmax_index
stats_quantile
stats_absdev
stats_absdev_m
stats_correlation
stats_covariance
stats_covariance_m
stats_kurtosis
stats_kurtosis_m_sd
stats_lag1_autocorrelation
stats_lag1_autocorrelation_m
stats_max
stats_max_index
stats_mean
stats_min
stats_min_index
stats_pvariance
stats_sd
stats_sd_m
stats_sd_with_fixed_mean
stats_skew
stats_skew_m_sd
stats_tss
stats_tss_m
stats_ttest
stats_variance
stats_variance_m
stats_variance_with_fixed_mean
stats_wabsdev
stats_wabsdev_m
stats_wkurtosis
stats_wkurtosis_m_sd
stats_wmean
stats_wsd
stats_wsd_m
stats_wsd_with_fixed_mean
stats_wskew
stats_wskew_m_sd
stats_wtss
stats_wtss_m
stats_wvariance
stats_wvariance_m
stats_wvariance_with_fixed_mean
);
%EXPORT_TAGS = ( 'all' => \@EXPORT_OK );

1;
