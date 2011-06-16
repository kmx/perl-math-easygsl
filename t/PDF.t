#!perl -T

use strict;
use warnings;

use Math::EasyGSL::PDF ':all';
use Test::Number::Delta tests => 35;
 
delta_within(pdf_gaussian(1.3, 0.6), 0.063588, 1e-4, 'pdf_gaussian');
delta_within(pdf_ugaussian(1.3), 0.171369, 1e-4, 'pdf_ugaussian');
delta_within(pdf_gaussian_tail(1.3, 1, 0.6), 1.330555, 1e-4, 'pdf_gaussian_tail');
delta_within(pdf_ugaussian_tail(1.3, 1), 1.080132, 1e-4, 'pdf_ugaussian_tail');
delta_within(pdf_bivariate_gaussian(1.8, 1.8, 2, 2, 0.9), 0.05960, 1e-4, 'pdf_bivariate_gaussian');;
delta_within(pdf_exponential(1.3, 2), 0.261023, 1e-4, 'pdf_exponential');
delta_within(pdf_laplace(1.3, 1), 0.13627, 1e-4, 'pdf_laplace');
delta_within(pdf_exppow(1.3, 1, 2), 0.10410, 1e-4, 'pdf_exppow');
delta_within(pdf_cauchy(1.3, 1), 0.11833, 1e-4, 'pdf_cauchy');
delta_within(pdf_rayleigh(1.3, 0.6), 0.34535, 1e-4, 'pdf_rayleigh');
delta_within(pdf_rayleigh_tail(1.3, 1, 0.6), 1.38498, 1e-4, 'pdf_rayleigh_tail');
delta_within(pdf_landau(1.3), 0.13251, 1e-4, 'pdf_landau');
delta_within(pdf_gamma(1.3, 1, 2), 0.26102, 1e-4, 'pdf_gamma');
delta_within(pdf_flat(1.3, 0, 2), 0.5, 1e-4, 'pdf_flat');
delta_within(pdf_lognormal(1.3, 0.4, 0.6), 0.49818, 1e-4, 'pdf_lognormal');
delta_within(pdf_chisq(1.3, 4), 0.16966, 1e-4, 'pdf_chisq');
delta_within(pdf_fdist(1.3, 4, 8), 0.32211, 1e-4, 'pdf_fdist');
delta_within(pdf_tdist(1.3, 4), 0.15538, 1e-4, 'pdf_tdist');
delta_within(pdf_beta(0.3, 11, 12), 0.90602, 1e-4, 'pdf_beta');
delta_within(pdf_logistic(1.3, 1), 0.16830, 1e-4, 'pdf_logistic');
delta_within(pdf_pareto(3, 1, 2), 0.22222, 1e-4, 'pdf_pareto');
delta_within(pdf_weibull(1.3, 1, 2), 0.47975, 1e-4, 'pdf_weibull');
delta_within(pdf_gumbel1(1.3, 1, 2), 0.31603, 1e-4, 'pdf_gumbel1');
delta_within(pdf_gumbel2(1.3, 1, 2), 0.25410, 1e-4, 'pdf_gumbel2');
delta_within(pdf_dirichlet(3, [1,4,5], [0.1,0.4,0.5]), 10.08, 1e-4, 'pdf_dirichlet');
delta_within(pdf_poisson(5, 4), 0.15629, 1e-4, 'pdf_poisson');
delta_within(pdf_bernoulli(5, 6), 0, 1e-4, 'pdf_bernoulli');
delta_within(pdf_binomial(25, 0.75, 33), 0.15943, 1e-4, 'pdf_binomial');
delta_within(pdf_multinomial(3, [1,2,3], [1,2,3]), 0.13889, 1e-4, 'pdf_multinomial');
delta_within(pdf_multinomial_ln(3, [1,2,3], [1,2,3]), -1.97408, 1e-4, 'pdf_multinomial_ln');
delta_within(pdf_negative_binomial(25, 0.65, 33), 0.02659, 1e-4, 'pdf_negative_binomial');
delta_within(pdf_pascal(5, 0.65, 3), 0.03029, 1e-4, 'pdf_pascal');
delta_within(pdf_geometric(5, 2.5), 12.65625, 1e-4, 'pdf_geometric');
delta_within(pdf_hypergeometric(5, 12, 18, 12), 0.29141, 1e-4, 'pdf_hypergeometric');
delta_within(pdf_logarithmic(2.5, 0.95), 0.15063, 1e-4, 'pdf_logarithmic');