use strict;
use warnings;
use ExtUtils::MakeMaker;

BEGIN {
  unless (`gsl-config --version` || `pkg-config --modversion gsl`) {
    warn "\n###FATAL### GSL library not found on your system\n";
    exit 0;
  }
}

WriteMakefile(
    NAME                => 'Math::EasyGSL',
    AUTHOR              => 'KMX <kmx@cpan.org>',
    VERSION_FROM        => 'lib/Math/EasyGSL.pm',
    ABSTRACT_FROM       => 'lib/Math/EasyGSL.pm',
    PREREQ_PM           => {
	'Exporter'      => 0,
        'Test::More'    => 0,
        'Carp'          => 0,
    },
    LIBS                => (`gsl-config --libs` || `pkg-config --libs gsl`),
    INC                 => (`gsl-config --cflags` || `pkg-config --cflags gsl`),
    NEEDS_LINKING       => 1,
    XS => {
      'lib/Math/EasyGSL/Statistics.xs' => 'lib/Math/EasyGSL/Statistics.c',
      'lib/Math/EasyGSL/Random.xs' => 'lib/Math/EasyGSL/Random.c',
      'lib/Math/EasyGSL/CDF.xs' => 'lib/Math/EasyGSL/CDF.c',      
      'lib/Math/EasyGSL/PDF.xs' => 'lib/Math/EasyGSL/PDF.c',
    },
    OBJECT => 'lib/Math/EasyGSL/Statistics.o lib/Math/EasyGSL/Random.o lib/Math/EasyGSL/CDF.o lib/Math/EasyGSL/PDF.o',
    META_MERGE => {
        resources => {
            license     => 'http://www.opensource.org/licenses/artistic-license-2.0.php',
            homepage    => 'http://github.com/kmx',
            bugtracker  => 'http://github.com/kmx',
            Repository  => 'http://github.com/kmx',
        }
    },
    LICENSE             => 'perl',
    MIN_PERL_VERSION    => 5.008,
    dist                => { COMPRESS => 'gzip -9f', SUFFIX => 'gz', },
    clean               => { FILES => 'Math-EasyGSL-*' },
);
