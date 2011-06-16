package Math::EasyGSL::PDF;

@ISA = qw/ DynaLoader /;
require DynaLoader;

bootstrap Math::EasyGSL::PDF;

use Exporter 'import'; # gives you Exporter's import() method directly
@EXPORT_OK = qw(
pdf_gaussian
pdf_ugaussian
pdf_gaussian_tail
pdf_ugaussian_tail
pdf_bivariate_gaussian
pdf_exponential
pdf_laplace
pdf_exppow
pdf_cauchy
pdf_rayleigh
pdf_rayleigh_tail
pdf_landau
pdf_gamma
pdf_flat
pdf_lognormal
pdf_chisq
pdf_fdist
pdf_tdist
pdf_beta
pdf_logistic
pdf_pareto
pdf_weibull
pdf_gumbel1
pdf_gumbel2
pdf_dirichlet
pdf_poisson
pdf_bernoulli
pdf_binomial
pdf_multinomial
pdf_multinomial_ln
pdf_negative_binomial
pdf_pascal
pdf_geometric
pdf_hypergeometric
pdf_logarithmic
);
%EXPORT_TAGS = ( 'all' => \@EXPORT_OK );

1;
