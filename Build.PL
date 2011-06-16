use strict;
use warnings;
use Module::Build;

BEGIN {
  unless (`gsl-config --version` || `pkg-config --modversion gsl`) {
    warn "\n###FATAL### GSL library not found on your system\n";
    exit 0;
  }
}

my $builder = Module::Build->new(
  module_name         => 'Math::EasyGSL',
  license             => 'perl',
  dist_author         => q{KMX <kmx@cpan.org>},
  dist_version_from   => 'lib/Math/EasyGSL.pm',
  build_requires      => { #need to have for running: ./Build (install|test)
      'Test::More'           => 0,
      'Test::Number::Delta'  => 0,
  },
  c_source            => 'src', # ppport.h
  create_readme => 1,
  extra_linker_flags  => `gsl-config --libs` || `pkg-config --libs gsl`,
  extra_compiler_flags=> `gsl-config --cflags` || `pkg-config --cflags gsl`,
  add_to_cleanup      => [ 'Math-EasyGSL-*' ],
  meta_merge          => {
    resources => {
      bugtracker => 'https://rt.cpan.org/NoAuth/Bugs.html?Dist=Math-EasyGSL',
      repository => 'https://github.com/kmx/perl-math-easygsl'
    },
  },
);

$builder->create_build_script();
