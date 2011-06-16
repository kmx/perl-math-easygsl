#!perl -T

use strict;
use warnings;

use Math::EasyGSL::Functions ':all';
use Test::Number::Delta tests => 2;
use Test::More;

delta_within(sf_bessel_Y0(1.23), 0.24638, 1e-4, 'sf_bessel_Y0');
delta_within([sf_bessel_Jn_array(5,8,4.44)], [0.18682,0.07942,0.02783,0.00833], 1e-4, 'sf_bessel_Jn_array');
