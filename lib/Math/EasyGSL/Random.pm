package Math::EasyGSL::Random;

@ISA = qw/ DynaLoader /;
require DynaLoader;

bootstrap Math::EasyGSL::Random;

sub new {
  my $class =shift;
  my %args = (
      'type' => undef,
      'env_setup' => 0,
      @_,
  );
  my $r = Math::EasyGSL::Random::_create($args{type}, $args{env_setup});
  my $self = { '!int!rnghandle' => $r };
  bless($self, $class);
  $self->seed($args{seed}) if $args{seed};
  return $self;
}

sub DESTROY {
  my $self = shift;
  my $r = Math::EasyGSL::Random::_destroy($self->{'!int!rnghandle'});
}

1;
