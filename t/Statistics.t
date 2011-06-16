#!perl -T

use strict;
use warnings;

use Math::EasyGSL::Random;
use Math::EasyGSL::Statistics ':all';

use Test::Number::Delta tests => 47;
use Test::More;

my $r = Math::EasyGSL::Random->new( seed=>5555 );

my @data  = map {$r->get_gaussian(0.5)} (1..100);
my @data1 = map {$r->get_gaussian(0.5)} (1..100);
my @data2 = map {$r->get_gaussian(0.5)} (1..100);
my @w = map {$r->get_gaussian(0.5)} (1..100);

delta_within(stats_median(\@data), -0.11479, 1e-4, 'stats_median');
my ($min, $max) = stats_minmax(\@data);
delta_within($min, -1.16677, 1e-4, 'stats_minmax (min)'); #xxx
delta_within($max, 1.42913, 1e-4, 'stats_minmax (max)'); #xxx
my ($min_i, $max_i) = stats_minmax_index(\@data);
is($min_i, 41, 'stats_minmax_index (min)');
is($max_i, 66, 'stats_minmax_index (max)');
delta_within(stats_quantile(\@data,0.18), -0.58971, 1e-4, 'stats_quantile');
delta_within(stats_absdev(\@data), 0.45729, 1e-4, 'stats_absdev');
delta_within(stats_absdev_m(\@data,0.5), 0.64121, 1e-4, 'stats_absdev_m');
delta_within(stats_correlation(\@data1,\@data2), -0.08671, 1e-4, 'stats_correlation');
delta_within(stats_covariance(\@data1,\@data2), -0.02070, 1e-4, 'stats_covariance');
delta_within(stats_covariance_m(\@data1,\@data2,0.4,0.5), 0.15284, 1e-4, 'stats_covariance_m');
delta_within(stats_kurtosis(\@data), -0.43873, 1e-4, 'stats_kurtosis');
delta_within(stats_kurtosis_m_sd(\@data,0.5,0.8), -1.14569, 1e-4, 'stats_kurtosis_m_sd');
delta_within(stats_lag1_autocorrelation(\@data), 0.05226, 1e-4, 'stats_lag1_autocorrelation');
delta_within(stats_lag1_autocorrelation_m(\@data,0.78), 0.67411, 1e-4, 'stats_lag1_autocorrelation_m');
delta_within(stats_max(\@data), 1.42913, 1e-4, 'stats_max');
is(stats_max_index(\@data), 66, 'stats_max_index');
delta_within(stats_mean(\@data), -0.02555, 1e-4, 'stats_mean');
delta_within(stats_min(\@data), -1.16677, 1e-4, 'stats_min');
is(stats_min_index(\@data), 41, 'stats_min_index');
delta_within(stats_pvariance(\@data1,\@data2), 0.23867, 1e-4, 'stats_pvariance');
delta_within(stats_sd(\@data), 0.56143, 1e-4, 'stats_sd');
delta_within(stats_sd_m(\@data,0.78), 0.98523, 1e-4, 'stats_sd_m');
delta_within(stats_sd_with_fixed_mean(\@data,0.78), 0.98029, 1e-4, 'stats_sd_with_fixed_mean');
delta_within(stats_skew(\@data), 0.23755, 1e-4, 'stats_skew');
delta_within(stats_skew_m_sd(\@data,0.78,0.2), -154.35169, 1e-4, 'stats_skew_m_sd');
delta_within(stats_tss(\@data), 31.20508, 1e-4, 'stats_tss');
delta_within(stats_tss_m(\@data,0.78), 96.09632, 1e-4, 'stats_tss_m');
delta_within(stats_ttest(\@data1,\@data2), 0.64646, 1e-4, 'stats_ttest');
delta_within(stats_variance(\@data), 0.31520, 1e-4, 'stats_variance');
delta_within(stats_variance_m(\@data,0.78), 0.97067, 1e-4, 'stats_variance_m');
delta_within(stats_variance_with_fixed_mean(\@data,0.78), 0.96096, 1e-4, 'stats_variance_with_fixed_mean');
delta_within(stats_wabsdev(\@w,\@data), 0.47396, 1e-4, 'stats_wabsdev');
delta_within(stats_wabsdev_m(\@w,\@data,0.88), 0.91048, 1e-4, 'stats_wabsdev_m');
delta_within(stats_wkurtosis(\@w,\@data), -0.73042, 1e-4, 'stats_wkurtosis');
delta_within(stats_wkurtosis_m_sd(\@w,\@data,0.88,0.21), 999.78438, 1e-4, 'stats_wkurtosis_m_sd');
delta_within(stats_wmean(\@w,\@data), 0.02592, 1e-4, 'stats_wmean');
delta_within(stats_wsd(\@w,\@data), 0.58972, 1e-4, 'stats_wsd');
delta_within(stats_wsd_m(\@w,\@data,0.88), 1.04805, 1e-4, 'stats_wsd_m');
delta_within(stats_wsd_with_fixed_mean(\@w,\@data,0.78), 0.95214, 1e-4, 'stats_wsd_with_fixed_mean');
delta_within(stats_wskew(\@w,\@data), 0.47878, 1e-4, 'stats_wskew');
delta_within(stats_wskew_m_sd(\@w,\@data,0.88,0.21), -150.16884, 1e-4, 'stats_wskew_m_sd');
delta_within(stats_wtss(\@w,\@data), 6.97796, 1e-4, 'stats_wtss');
delta_within(stats_wtss_m(\@w,\@data,0.88), 22.03944, 1e-4, 'stats_wtss_m');
delta_within(stats_wvariance(\@w,\@data), 0.34777, 1e-4, 'stats_wvariance');
delta_within(stats_wvariance_m(\@w,\@data,0.88), 1.09841, 1e-4, 'stats_wvariance_m');
delta_within(stats_wvariance_with_fixed_mean(\@w,\@data,0.78), 0.90658, 1e-4, 'stats_wvariance_with_fixed_mean');
