package Math::EasyGSL::CDF;

@ISA = qw/ DynaLoader /;
require DynaLoader;

use Exporter 'import'; # gives you Exporter's import() method directly
@EXPORT_OK = qw(
cdf_beta_P
cdf_beta_Pinv
cdf_beta_Q
cdf_beta_Qinv
cdf_binomial_P
cdf_binomial_Q
cdf_cauchy_P
cdf_cauchy_Pinv
cdf_cauchy_Q
cdf_cauchy_Qinv
cdf_chisq_P
cdf_chisq_Pinv
cdf_chisq_Q
cdf_chisq_Qinv
cdf_exponential_P
cdf_exponential_Pinv
cdf_exponential_Q
cdf_exponential_Qinv
cdf_exppow_P
cdf_exppow_Q
cdf_fdist_P
cdf_fdist_Pinv
cdf_fdist_Q
cdf_fdist_Qinv
cdf_flat_P
cdf_flat_Pinv
cdf_flat_Q
cdf_flat_Qinv
cdf_gamma_P
cdf_gamma_Pinv
cdf_gamma_Q
cdf_gamma_Qinv
cdf_gaussian_P
cdf_gaussian_Pinv
cdf_gaussian_Q
cdf_gaussian_Qinv
cdf_geometric_P
cdf_geometric_Q
cdf_gumbel1_P
cdf_gumbel1_Pinv
cdf_gumbel1_Q
cdf_gumbel1_Qinv
cdf_gumbel2_P
cdf_gumbel2_Pinv
cdf_gumbel2_Q
cdf_gumbel2_Qinv
cdf_hypergeometric_P
cdf_hypergeometric_Q
cdf_laplace_P
cdf_laplace_Pinv
cdf_laplace_Q
cdf_laplace_Qinv
cdf_logistic_P
cdf_logistic_Pinv
cdf_logistic_Q
cdf_logistic_Qinv
cdf_lognormal_P
cdf_lognormal_Pinv
cdf_lognormal_Q
cdf_lognormal_Qinv
cdf_negative_binomial_P
cdf_negative_binomial_Q
cdf_pareto_P
cdf_pareto_Pinv
cdf_pareto_Q
cdf_pareto_Qinv
cdf_pascal_P
cdf_pascal_Q
cdf_poisson_P
cdf_poisson_Q
cdf_rayleigh_P
cdf_rayleigh_Pinv
cdf_rayleigh_Q
cdf_rayleigh_Qinv
cdf_tdist_P
cdf_tdist_Pinv
cdf_tdist_Q
cdf_tdist_Qinv
cdf_ugaussian_P
cdf_ugaussian_Pinv
cdf_ugaussian_Q
cdf_ugaussian_Qinv
cdf_weibull_P
cdf_weibull_Pinv
cdf_weibull_Q
cdf_weibull_Qinv
);
%EXPORT_TAGS = ( 'all' => \@EXPORT_OK );

bootstrap Math::EasyGSL::CDF;

1;
