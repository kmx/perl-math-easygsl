#!perl -T

use strict;
use warnings;

use Math::EasyGSL ':all';
use Test::More tests => 3;

like(GSL_VERSION, qr/\d+\.\d+/, 'GSL_VERSION');
cmp_ok(GSL_MAJOR_VERSION, '>=', 1, 'GSL_MAJOR_VERSION');
cmp_ok(GSL_MINOR_VERSION, '>=', 0, 'GSL_MINOR_VERSION');
