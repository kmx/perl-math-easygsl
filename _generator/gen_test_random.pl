use strict;
use warnings;

my @l = (
'get_gaussian(0.6)',
'get_ugaussian()',
'get_gaussian_tail(1, 0.6)',
'get_ugaussian_tail(1)',
'get_bivariate_gaussian(2, 2, 0.9)',
'get_exponential(2)',
'get_laplace(1)',
'get_exppow(1, 2)',
'get_cauchy(1)',
'get_rayleigh(0.6)',
'get_rayleigh_tail(1, 0.6)',
'get_landau()',
'get_gamma(1, 2)',
'get_flat(0, 2)',
'get_lognormal(0.4, 0.6)',
'get_chisq(4)',
'get_fdist(4, 8)',
'get_tdist(4)',
'get_beta(11, 12)',
'get_logistic(1)',
'get_pareto(1, 2)',
'get_weibull(1, 2)',
'get_gumbel1(1, 2)',
'get_gumbel2(1, 2)',
#'get_dirichlet([1,4,5], [0.1,0.4,0.5])',
'get_poisson(4)',
'get_bernoulli(6)',
'get_binomial(0.75, 33)',
#'get_multinomial([1,2,3], [1,2,3])',
#'get_multinomial_ln([1,2,3], [1,2,3])',
'get_negative_binomial(0.65, 33)',
'get_pascal(0.65, 3)',
'get_geometric(2.5)',
'get_hypergeometric(12, 18, 12)',
'get_logarithmic(0.95)',
);

use Math::EasyGSL::Random;
my $r = Math::EasyGSL::Random->new( seed=>5555 );

my $output = <<'MARKER';
use strict;
use warnings;
use Math::EasyGSL::Random;
my $r = Math::EasyGSL::Random->new( seed=>5555 );
MARKER

for (@l) {
  if (/^([^\(]+)(.*)/) {
    #print "1=$1 2=$2\n";
    $output .= sprintf q!print 'delta_within([map {$r->%s%s} (1..5)], ['.join(',',map {sprintf("%%.5f",$_)} (map {$r->%s%s} (1..5)))."], 1e-4, '%s');\n";!, $1,$2,$1,$2,$1;
    $output .= "\n";
  }
  else {
    die "Invalid item '$_'\n";
  }
}

eval $output;
warn $@ if $@;