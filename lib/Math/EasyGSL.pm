package Math::EasyGSL;

=head1 NAME

Math::EasyGSL - The great new Math::EasyGSL!

=head1 VERSION

Version 0.01

=cut

our $VERSION = "0.01";

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

sub xxx {
}

=head1 SYNOPSIS

Perhaps a little code snippet.

    use Math::EasyGSL;

    my $foo = Math::EasyGSL->new();
    ...

=head1 EXPORT

A list of functions that can be exported.  You can delete this section
if you don't export anything, such as for a purely object-oriented module.

=head1 SUBROUTINES/METHODS

=head1 LICENSE AND COPYRIGHT

Copyright 2011 KMX.

This program is free software; you can redistribute it and/or modify it
under the terms of either: the GNU General Public License as published
by the Free Software Foundation; or the Artistic License.

See http://dev.perl.org/licenses/ for more information.

=cut

1; # End of Math::EasyGSL