#!perl -T

use strict;
use warnings;

use Math::EasyGSL::CDF ':all';
use Test::Number::Delta tests => 86;
 
delta_within(cdf_gaussian_P(0.5, 0.6), 0.797672, 1e-4, 'cdf_gaussian_P');
delta_within(cdf_gaussian_Q(0.5, 0.6), 0.202328, 1e-4, 'cdf_gaussian_Q');
delta_within(cdf_gaussian_Pinv(0.7, 0.6), 0.31464, 1e-4, 'cdf_gaussian_Pinv');
delta_within(cdf_gaussian_Qinv(0.8, 0.6), -0.50497, 1e-4, 'cdf_gaussian_Qinv');
delta_within(cdf_ugaussian_P(0.5), 0.69146, 1e-4, 'cdf_ugaussian_P');
delta_within(cdf_ugaussian_Q(0.5), 0.30854, 1e-4, 'cdf_ugaussian_Q');
delta_within(cdf_ugaussian_Pinv(0.7), 0.52440, 1e-4, 'cdf_ugaussian_Pinv');
delta_within(cdf_ugaussian_Qinv(0.8), -0.84162, 1e-4, 'cdf_ugaussian_Qinv');
delta_within(cdf_exponential_Pinv(0.7, 0.25), 0.30099, 1e-4, 'cdf_exponential_Pinv');
delta_within(cdf_exponential_Q(0.5, 0.25), 0.13534, 1e-4, 'cdf_exponential_Q');
delta_within(cdf_exponential_P(0.5, 0.25), 0.86466, 1e-4, 'cdf_exponential_P');
delta_within(cdf_exponential_Qinv(0.8, 0.25), 0.05579, 1e-4, 'cdf_exponential_Qinv');
delta_within(cdf_laplace_P(0.5, 100), 0.50249, 1e-4, 'cdf_laplace_P');
delta_within(cdf_laplace_Q(0.5, 100), 0.49751, 1e-4, 'cdf_laplace_Q');
delta_within(cdf_laplace_Pinv(0.7, 100), 51.08256, 1e-4, 'cdf_laplace_Pinv');
delta_within(cdf_laplace_Qinv(0.8, 100), -91.62907, 1e-4, 'cdf_laplace_Qinv');
delta_within(cdf_exppow_P(0.5, 10, 20), 0.52568, 1e-4, 'cdf_exppow_P');
delta_within(cdf_exppow_Q(0.5, 10, 20), 0.47432, 1e-4, 'cdf_exppow_Q');
delta_within(cdf_cauchy_P(0.5, 10), 0.51590, 1e-4, 'cdf_cauchy_P');
delta_within(cdf_cauchy_Q(0.5, 10), 0.48410, 1e-4, 'cdf_cauchy_Q');
delta_within(cdf_cauchy_Pinv(0.7, 10), 7.26543, 1e-4, 'cdf_cauchy_Pinv');
delta_within(cdf_cauchy_Qinv(0.8, 10), -13.76382, 1e-4, 'cdf_cauchy_Qinv');
delta_within(cdf_rayleigh_P(0.5, 0.6), 0.29335, 1e-4, 'cdf_rayleigh_P');
delta_within(cdf_rayleigh_Q(0.5, 0.6), 0.70665, 1e-4, 'cdf_rayleigh_Q');
delta_within(cdf_rayleigh_Pinv(0.7, 0.6), 0.93105, 1e-4, 'cdf_rayleigh_Pinv');
delta_within(cdf_rayleigh_Qinv(0.8, 0.6), 0.40083, 1e-4, 'cdf_rayleigh_Qinv');
delta_within(cdf_gamma_P(0.5, 1, 2), 0.22120, 1e-4, 'cdf_gamma_P');
delta_within(cdf_gamma_Q(0.5, 1, 2), 0.77880, 1e-4, 'cdf_gamma_Q');
delta_within(cdf_gamma_Pinv(0.7, 1, 2), 2.40795, 1e-4, 'cdf_gamma_Pinv');
delta_within(cdf_gamma_Qinv(0.8, 1, 2), 0.44629, 1e-4, 'cdf_gamma_Qinv');
delta_within(cdf_flat_P(0.5, 0, 2), 0.25, 1e-4, 'cdf_flat_P');
delta_within(cdf_flat_Q(0.5, 0, 2), 0.75, 1e-4, 'cdf_flat_Q');
delta_within(cdf_flat_Pinv(0.7, 1, 2), 1.7, 1e-4, 'cdf_flat_Pinv');
delta_within(cdf_flat_Qinv(0.8, 1, 2), 1.2, 1e-4, 'cdf_flat_Qinv');
delta_within(cdf_lognormal_P(0.5, 0.3, 0.6), 0.04894, 1e-4, 'cdf_lognormal_P');
delta_within(cdf_lognormal_Q(0.5, 0.3, 0.6), 0.95106, 1e-4, 'cdf_lognormal_Q');
delta_within(cdf_lognormal_Pinv(0.7, 0.3, 0.6), 1.84899, 1e-4, 'cdf_lognormal_Pinv');
delta_within(cdf_lognormal_Qinv(0.8, 0.3, 0.6), 0.81467, 1e-4, 'cdf_lognormal_Qinv');
delta_within(cdf_chisq_P(0.5, 0.9), 0.56140, 1e-4, 'cdf_chisq_P');
delta_within(cdf_chisq_Q(0.5, 0.9), 0.43860, 1e-4, 'cdf_chisq_Q');
delta_within(cdf_chisq_Pinv(0.7, 0.9), 0.93483, 1e-4, 'cdf_chisq_Pinv');
delta_within(cdf_chisq_Qinv(0.8, 0.9), 0.04336, 1e-4, 'cdf_chisq_Qinv');
delta_within(cdf_fdist_P(0.5, 0.5, 0.8), 0.49658, 1e-4, 'cdf_fdist_P');
delta_within(cdf_fdist_Q(0.5, 0.5, 0.8), 0.50342, 1e-4, 'cdf_fdist_Q');
delta_within(cdf_fdist_Pinv(0.7, 1.5, 1.8), 2.36542, 1e-4, 'cdf_fdist_Pinv');
delta_within(cdf_fdist_Qinv(0.8, 1.5, 1.8), 0.17989, 1e-4, 'cdf_fdist_Qinv');
delta_within(cdf_tdist_P(0.5, 0.7), 0.63480, 1e-4, 'cdf_tdist_P');
delta_within(cdf_tdist_Q(0.5, 0.7), 0.36520, 1e-4, 'cdf_tdist_Q');
delta_within(cdf_tdist_Pinv(0.7, 2), 0.61721, 1e-4, 'cdf_tdist_Pinv');
delta_within(cdf_tdist_Qinv(0.8, 2), -1.06066, 1e-4, 'cdf_tdist_Qinv');
delta_within(cdf_beta_P(0.5, 1, 2), 0.75, 1e-4, 'cdf_beta_P');
delta_within(cdf_beta_Q(0.5, 1, 2), 0.25, 1e-4, 'cdf_beta_Q');
delta_within(cdf_beta_Pinv(0.7, 1, 2), 0.45228, 1e-4, 'cdf_beta_Pinv');
delta_within(cdf_beta_Qinv(0.8, 1, 2), 0.10557, 1e-4, 'cdf_beta_Qinv');
delta_within(cdf_logistic_P(0.5, 2), 0.56218, 1e-4, 'cdf_logistic_P');
delta_within(cdf_logistic_Q(0.5, 2), 0.43782, 1e-4, 'cdf_logistic_Q');
delta_within(cdf_logistic_Pinv(0.7, 2), 1.69460, 1e-4, 'cdf_logistic_Pinv');
delta_within(cdf_logistic_Qinv(0.8, 2), -2.77259, 1e-4, 'cdf_logistic_Qinv');
delta_within(cdf_pareto_P(4.9, 1, 2), 0.59184, 1e-4, 'cdf_pareto_P');
delta_within(cdf_pareto_Q(4.9, 1, 2), 0.40816, 1e-4, 'cdf_pareto_Q');
delta_within(cdf_pareto_Pinv(0.3, 11, 22), 22.72504, 1e-4, 'cdf_pareto_Pinv');
delta_within(cdf_pareto_Qinv(0.3, 11, 22), 24.54467, 1e-4, 'cdf_pareto_Qinv');
delta_within(cdf_weibull_P(0.5, 1, 2), 0.22120, 1e-4, 'cdf_weibull_P');
delta_within(cdf_weibull_Q(0.5, 1, 2), 0.77880, 1e-4, 'cdf_weibull_Q');
delta_within(cdf_weibull_Pinv(0.7, 1, 2), 1.09726, 1e-4, 'cdf_weibull_Pinv');
delta_within(cdf_weibull_Qinv(0.8, 1, 2), 0.47238, 1e-4, 'cdf_weibull_Qinv');
delta_within(cdf_gumbel1_P(0.5, 1, 2), 0.29729, 1e-4, 'cdf_gumbel1_P');
delta_within(cdf_gumbel1_Q(0.5, 1, 2), 0.70271, 1e-4, 'cdf_gumbel1_Q');
delta_within(cdf_gumbel1_Pinv(0.7, 1, 2), 1.72408, 1e-4, 'cdf_gumbel1_Pinv');
delta_within(cdf_gumbel1_Qinv(0.8, 1, 2), 0.21726, 1e-4, 'cdf_gumbel1_Qinv');
delta_within(cdf_gumbel2_P(0.5, 1, 2), 0.01832, 1e-4, 'cdf_gumbel2_P');
delta_within(cdf_gumbel2_Q(0.5, 1, 2), 0.98168, 1e-4, 'cdf_gumbel2_Q');
delta_within(cdf_gumbel2_Pinv(0.7, 1, 2), 5.60735, 1e-4, 'cdf_gumbel2_Pinv');
delta_within(cdf_gumbel2_Qinv(0.8, 1, 2), 1.24267, 1e-4, 'cdf_gumbel2_Qinv');
delta_within(cdf_poisson_P(2, 0.5), 0.98561, 1e-4, 'cdf_poisson_P');
delta_within(cdf_poisson_Q(2, 0.5), 0.01439, 1e-4, 'cdf_poisson_Q');
delta_within(cdf_binomial_P(11, 0.3, 33), 0.73338, 1e-4, 'cdf_binomial_P');
delta_within(cdf_binomial_Q(11, 0.3, 33), 0.26662, 1e-4, 'cdf_binomial_Q');
delta_within(cdf_negative_binomial_P(11, 0.3, 3), 0.83916, 1e-4, 'cdf_negative_binomial_P');
delta_within(cdf_negative_binomial_Q(11, 0.3, 3), 0.16084, 1e-4, 'cdf_negative_binomial_Q');
delta_within(cdf_pascal_P(11, 0.3, 3), 0.83916, 1e-4, 'cdf_pascal_P');
delta_within(cdf_pascal_Q(11, 0.3, 3), 0.16084, 1e-4, 'cdf_pascal_Q');
delta_within(cdf_geometric_P(11, 0.3), 0.98023, 1e-4, 'cdf_geometric_P');
delta_within(cdf_geometric_Q(11, 0.3), 0.01977, 1e-4, 'cdf_geometric_Q');
delta_within(cdf_hypergeometric_P(2, 50, 2, 3), 0.11312, 1e-4, 'cdf_hypergeometric_P');
delta_within(cdf_hypergeometric_Q(2, 50, 2, 3), 0.88688, 1e-4, 'cdf_hypergeometric_Q');
