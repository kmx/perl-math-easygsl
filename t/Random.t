#!perl -T

use strict;
use warnings;

use Math::EasyGSL::Random;
use Test::Number::Delta tests => 39;
use Test::More;
 
my $r = Math::EasyGSL::Random->new( seed=>5555 );

delta_within([map {$r->get_gaussian(0.6)} (1..5)], [-0.74561,0.52873,0.48301,-0.40675,-0.19424], 1e-4, 'get_gaussian');
delta_within([map {$r->get_ugaussian()} (1..5)], [0.12691,0.39379,-0.46870,-0.12130,1.74730], 1e-4, 'get_ugaussian');
delta_within([map {$r->get_gaussian_tail(1, 0.6)} (1..5)], [1.24115,1.39139,1.03054,1.25289,1.11175], 1e-4, 'get_gaussian_tail');
delta_within([map {$r->get_ugaussian_tail(1)} (1..5)], [1.95449,1.61052,1.38371,1.35216,2.01014], 1e-4, 'get_ugaussian_tail');
delta_within([map {$r->get_bivariate_gaussian(2, 2, 0.9)} (1..5)], [0.59960,0.33302,-0.50370,-1.01061,2.20250,1.41779,-0.51444,-0.37362,1.64632,1.31046], 1e-4, 'get_bivariate_gaussian');
delta_within([map {$r->get_exponential(2)} (1..5)], [3.32521,0.42790,0.13040,4.55218,1.14462], 1e-4, 'get_exponential');
delta_within([map {$r->get_laplace(1)} (1..5)], [2.33438,1.36069,0.57335,0.64465,0.71160], 1e-4, 'get_laplace');
delta_within([map {$r->get_exppow(1, 2)} (1..5)], [0.33953,-0.25436,-1.65006,-0.68252,0.86250], 1e-4, 'get_exppow');
delta_within([map {$r->get_cauchy(1)} (1..5)], [-1.67820,-4.93554,-0.59154,-0.31480,2.02054], 1e-4, 'get_cauchy');
delta_within([map {$r->get_rayleigh(0.6)} (1..5)], [0.17739,0.66113,0.86694,0.57002,0.62539], 1e-4, 'get_rayleigh');
delta_within([map {$r->get_rayleigh_tail(1, 0.6)} (1..5)], [1.03959,1.48064,1.39000,1.79004,1.37949], 1e-4, 'get_rayleigh_tail');
delta_within([map {$r->get_landau()} (1..5)], [1.10767,-0.06387,0.38527,2.11301,0.18270], 1e-4, 'get_landau');
delta_within([map {$r->get_gamma(1, 2)} (1..5)], [0.49401,0.64388,0.26392,0.38014,1.22863], 1e-4, 'get_gamma');
delta_within([map {$r->get_flat(0, 2)} (1..5)], [0.54023,0.61952,0.58681,0.96064,0.25827], 1e-4, 'get_flat');
delta_within([map {$r->get_lognormal(0.4, 0.6)} (1..5)], [1.31388,1.79322,0.60784,0.61469,0.91616], 1e-4, 'get_lognormal');
delta_within([map {$r->get_chisq(4)} (1..5)], [1.97381,4.16687,1.06232,4.12531,6.50126], 1e-4, 'get_chisq');
delta_within([map {$r->get_fdist(4, 8)} (1..5)], [0.65367,0.43243,0.60507,1.27071,0.54377], 1e-4, 'get_fdist');
delta_within([map {$r->get_tdist(4)} (1..5)], [-2.96538,0.35844,0.28280,0.48989,-3.67999], 1e-4, 'get_tdist');
delta_within([map {$r->get_beta(11, 12)} (1..5)], [0.60249,0.48429,0.22454,0.56026,0.48175], 1e-4, 'get_beta');
delta_within([map {$r->get_logistic(1)} (1..5)], [-0.25211,1.62553,0.66153,1.06243,2.92398], 1e-4, 'get_logistic');
delta_within([map {$r->get_pareto(1, 2)} (1..5)], [38.20978,15.28055,3.93719,4.82515,3.58857], 1e-4, 'get_pareto');
delta_within([map {$r->get_weibull(1, 2)} (1..5)], [0.35316,0.28707,1.24270,0.89751,0.67818], 1e-4, 'get_weibull');
delta_within([map {$r->get_gumbel1(1, 2)} (1..5)], [0.63851,1.31749,1.69514,2.09542,0.51708], 1e-4, 'get_gumbel1');
delta_within([map {$r->get_gumbel2(1, 2)} (1..5)], [2.16792,1.62736,1.20844,12.13493,1.08872], 1e-4, 'get_gumbel2');
delta_within([map {$r->get_poisson(4)} (1..5)], [5,2,3,6,2], 1e-4, 'get_poisson');
delta_within([map {$r->get_bernoulli(6)} (1..5)], [1,1,1,1,1], 1e-4, 'get_bernoulli');
delta_within([map {$r->get_binomial(0.75, 33)} (1..5)], [23,27,22,24,24], 1e-4, 'get_binomial');
delta_within([map {$r->get_negative_binomial(0.65, 33)} (1..5)], [9,16,16,13,24], 1e-4, 'get_negative_binomial');
delta_within([map {$r->get_pascal(0.65, 3)} (1..5)], [0,3,2,1,1], 1e-4, 'get_pascal');
delta_within([map {$r->get_geometric(2.5)} (1..5)], [0,0,0,0,0], 1e-4, 'get_geometric');
delta_within([map {$r->get_hypergeometric(12, 18, 12)} (1..5)], [6,4,3,4,6], 1e-4, 'get_hypergeometric');
delta_within([map {$r->get_logarithmic(0.95)} (1..5)], [2,3,1,1,1], 1e-4, 'get_logarithmic');

delta_within( $r->get_dirichlet([1,2]), [0.33077,0.66923], 1e-4, 'get_dirichlet');
delta_within( $r->get_multinomial(3, [0.15,0.15,0.3,0.2,0.2]), [0,0,2,1,0], 1e-4, 'get_multinomial');
is_deeply( $r->get_poisson_array(3, 2), [3,1,4], 'get_poisson_array');
is($r->get_discrete([0.05,0.05,0.6,0.1,0.2]), 2, 'get_discrete');

is_deeply( $r->shuffle([44,55,66,77,88,99]), [44,88,55,99,77,66], 'shuffle');
is_deeply( $r->choose(3,[44,55,66,77,88,99]), [44,55,77], 'choose');
is_deeply( $r->sample(5,[7,8]), [7,7,8,8,8], 'sample');
