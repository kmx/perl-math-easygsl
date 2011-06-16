package Math::EasyGSL;

=head1 NAME

Math::EasyGSL - Perl bindings to GSL (GNU Scientific Library)

=head1 VERSION

Version 0.001

=cut

# following recommendation from http://www.dagolden.com/index.php/369/version-numbers-should-be-boring/
our $VERSION = "0.001";
$VERSION = eval $VERSION;

@ISA = qw/ DynaLoader /;
require DynaLoader;

bootstrap Math::EasyGSL;

use Exporter 'import'; # gives you Exporter's import() method directly
@EXPORT_OK = qw(
GSL_VERSION
GSL_MINOR_VERSION
GSL_MAJOR_VERSION
);
%EXPORT_TAGS = ( 'all' => \@EXPORT_OK );

=head1 SYNOPSIS

 # the main module Math::EasyGSL has just functions providing simple info
 use Math::EasyGSL ':all';
 print "Welcome to Math::EasyGSL version=" . GSL_VERSION;
 
 # the real math related function are in accompanying modules
 use Math::EasyGSL::Statistics ':all';
 my @data = (11, 12, 13, 14, 15);
 print "Mean=" . stats_mean(\@data);

This module B<is not intended> as a replacement for L<Math::GSL>.

Math::EasyGSL is not covering 100% of all GSL functions. It binds into
perl just part of GSL library (I am adding more functions on "if-needed"
basis - if you miss something let me know).

On the other hand the interface of L<Math::EasyGSL> is more perl-friendly than
L<Math::GSL> - it is not 1:1 translation of the original C functions.
 
=head1 EXPORT

By default there are no functions exported. You can import into your
program all available functions by calling:

 use Math::EasyGSL ':all';

Or you can import just selected functions:

 use Math::EasyGSL::PDF qw(GSL_VERSION);

There are no other import tags.

=head1 SUBROUTINES/METHODS

=head2 Version related functions

=head3 GSL_VERSION

Returns a string with version of underlaying GSL library like: C<1.14>

=head3 GSL_MAJOR_VERSION

Returns major version number like: C<1> (for GSL v1.14)

=head3 GSL_MINOR_VERSION

Returns version number like: C<14> (for GSL v1.14)

=head1 LICENSE AND COPYRIGHT

Copyright 2011 KMX.

This program is free software; you can redistribute it and/or modify it
under the terms of either: the GNU General Public License as published
by the Free Software Foundation; or the Artistic License.

See http://dev.perl.org/licenses/ for more information.

=cut

1;