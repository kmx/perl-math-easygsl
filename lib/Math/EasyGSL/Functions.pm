package Math::EasyGSL::Functions;

@ISA = qw/ DynaLoader /;
require DynaLoader;

bootstrap Math::EasyGSL::Functions;

use Exporter 'import'; # gives you Exporter's import() method directly
@EXPORT_OK = qw(
        sf_bessel_Y0
        sf_bessel_Jn_array
);
%EXPORT_TAGS = ( 'all' => \@EXPORT_OK );

1;
