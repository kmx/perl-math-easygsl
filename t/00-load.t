#!perl -T

use Test::More tests => 5;

BEGIN {
    use_ok( 'Math::EasyGSL' )             || print "Failed to use Math::EasyGSL!";
    use_ok( 'Math::EasyGSL::Random' )     || print "Failed to use Math::EasyGSL::Random!";
    use_ok( 'Math::EasyGSL::PDF' )        || print "Failed to use Math::EasyGSL::PDF!";
    use_ok( 'Math::EasyGSL::CDF' )        || print "Failed to use Math::EasyGSL::CDF!";
    use_ok( 'Math::EasyGSL::Statistics' ) || print "Failed to use Math::EasyGSL::Statistics!";
}

diag( "Testing Math::EasyGSL $Math::EasyGSL::VERSION, Perl $], $^X" );
