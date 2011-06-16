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
    build_requires => {
        'Test::More' => 0,
    },
    c_source            => 'src', # ppport.h
    extra_linker_flags  => `gsl-config --libs` || `pkg-config --libs gsl`,
    extra_compiler_flags=> `gsl-config --cflags` || `pkg-config --cflags gsl`,  
    add_to_cleanup      => [ 'Math-EasyGSL-*' ],
);

$builder->create_build_script();
