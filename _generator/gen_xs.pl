use strict;
use warnings;

use Template;
use File::Spec::Functions qw/rel2abs/;
use Data::Dumper;

sub proc_list {
  my ($intro, $list, $qr, $template, $outxs) = @_; 
  my $allinone = {};
  my @allfnc;

  warn "### proc_list('$list')\n";
  if (-f $list) {
    open DAT, "<", $list;
    while (<DAT>) {
      s/[\r\n]*$//;
      next if /^\s*$/;
      next if /^\/\//;
      #warn "re=$qr tx='$_'\n";
      next unless /$qr/;
      
      if (/^.*?#([^#]*)#(.*?)\s+(gsl_[a-z]+_[^\s\(]+)\s*\(([^\)]*)\)/) {
        #warn "1='$1' 2='$2' 3='$3' 4='$4'\n";	
	my $orig = "$2 $3($4);";	
	my ($xs_name, $gsl_rv, $gsl_name) = ($1, $2, $3);
	warn "func=$xs_name\n";
        my @pars = split ',', $4;
        s/^ *// for @pars;
        s/ *$// for @pars;
	
	my $svfunc;
	$svfunc = 'newSVnv' if $gsl_rv eq 'double';	
	$svfunc = 'newSViv' if $gsl_rv =~ /(size_t|unsigned int|int|unsigned long|long)/;
	$svfunc = 'XXX_void_XXX' if $gsl_rv eq 'void';
	die "ERROR: UNKNOWN svfunc for '$gsl_rv'\n" unless $svfunc;
	die "ERROR: unexpected gsl_rv=$gsl_rv" unless $gsl_rv;
	
	my $r2val;
	$r2val = 'r2double' if $orig =~ /gsl_sf_result/ && $gsl_name =~ /_e$/;
	$r2val = 'e10r2double' if $orig =~ /gsl_sf_result_e10/ && $gsl_name =~ /_e$/;
	my $result_array_size;
	my $result_array_size_msg;
	
	my @xs_params;
	my @gsl_params;
	my @xs_partypes;
	my @init_decl;
	my @init_arr;
	my @init_arr_size;
	my @init_arr_size_if;
	my $n_first;
	my @n_other;
	my @xs_args;

	for (@pars) {
	  #warn "p='$_'\n";
          if (/^(.*?)([a-zA-Z0-9_]*)([\[\]]*)?$/) {
            my ($n, $t) = ($2, $1);
	    my $arr = $3 ? 1 : 0; # [] in the end
	    $t =~ s/\s*$//;
	    $t =~ s/^\s*//;
	    my $ptr = ($t =~ /\*\s*/) ? 1 : 0;
	    $t =~ s/^const\s*//;	    
	    $n =~ s/\[\]\s*$//;
	    #warn "'$n'~'$t' arr=$arr ptr=$ptr\n";
	    if ($n eq 'r' && $t eq 'gsl_rng *') {
	      push @xs_args, { name=>'self', type=>'SV *', decl=>'SV * self' };
	      push @xs_params, 'self';
	      push @xs_partypes, 'SV * self';
	      push @gsl_params, 'ref2rng(self)';
	    }
	    elsif ($t eq 'int' && $n eq 'nmin' && $gsl_name =~ /_array/ && !defined($result_array_size_msg)) {
	      $result_array_size = 'nmax-nmin+1';
	      $result_array_size_msg = "invalid parameter 'nmax' and/or 'nmin'";
	      push @xs_args, { name=>$n, type=>$t, decl=>"$t $n" };
	      push @xs_params, $n;
	      push @xs_partypes, "$t $n";
	      push @gsl_params, $n;
	    }
	    elsif ($t eq 'int' && $n =~ /^(nmax|lmax)$/ && $gsl_name =~ /_array/ && !defined($result_array_size_msg)) {
	      $result_array_size = "$n+1";
	      $result_array_size_msg = "invalid parameter '$n'";
	      push @xs_args, { name=>$n, type=>$t, decl=>"$t $n" };
	      push @xs_params, $n;
	      push @xs_partypes, "$t $n";
	      push @gsl_params, $n;
	    }
	    elsif ($t eq 'int' && $n eq 'lmax' && $gsl_name =~ /_array/ && !defined($result_array_size_msg)) {
	      $result_array_size = 'nmax+1';
	      $result_array_size_msg = "invalid parameter 'nmax'";
	      push @xs_args, { name=>$n, type=>$t, decl=>"$t $n" };
	      push @xs_params, $n;
	      push @xs_partypes, "$t $n";
	      push @gsl_params, $n;
	    }
	    elsif ($t eq 'double *' && $n =~ /^(jl_x|result)_array/) {
	      push @gsl_params, 'result_array_';
	      push @init_decl, "double * result_array_";
	      push @init_decl, "int i, result_array_size";
	    }
	    elsif ($t eq 'gsl_sf_result *') {
	      push @gsl_params, "\&$n\_";	    
	      push @init_decl, "gsl_sf_result $n\_";
	      #push @init_arr, { name=>"$n\_", type=>'gsl_sf_result', origname=>$n, sv2val=>'xxxxx' };
	    }
	    elsif ($t eq 'gsl_sf_result_e10 *') {
	      push @gsl_params, "\&$n\_";	    
	      push @init_decl, "gsl_sf_result_e10 $n\_";
	      #push @init_arr, { name=>"$n\_", type=>'gsl_sf_result_e10', origname=>$n, sv2val=>'xxxxx' };
	    }
	    elsif ($arr) {
	      push @xs_args, { name=>$n, type=>'SV *', decl=>"SV * $n" };
	      push @xs_params, $n;
	      push @xs_partypes, "SV * $n";
	      push @gsl_params, "$n\_";	    
	      push @init_decl, "$t * $n\_";
	      push @init_decl, "size_t $n\_size";	    
	      push @init_arr, { name=>"$n\_", type=>$t, origname=>$n, sv2val=>($t=~/double/?'SvNV':'SvIV') };
	      push @init_arr_size, "$n\_size = arr_ref_size($n, \"Warning: $xs_name() - invalid '$n' argument\")";
	      push @init_arr_size_if, "$n\_size";
	      if ($n_first) {
	        push @n_other, "$n\_size";
	      }
	      else {
	        $n_first = "$n\_size";
	      }
	    }
	    elsif ($t eq 'gsl_mode_t' && $n eq 'mode') {
	      push @xs_args, { name=>$n, type=>'unsigned int', decl=>"unsigned int $n" };
	      push @xs_params, $n;
	      push @xs_partypes, "unsigned int $n";
	      push @gsl_params, $n;
	    }
	    else {	      
	      push @xs_args, { name=>$n, type=>$t, decl=>"$t $n" };
	      push @xs_params, $n;
	      push @xs_partypes, "$t $n";
	      push @gsl_params, $n;
	    }
          }
        }

	my @xs_args_censored;
	if (scalar(@init_arr)) {
	  @gsl_params = map { (/^w?stride([1-9]*)$/)?1:$_ } @gsl_params;
	  @gsl_params = map { (/^n([1-9]*)$/)?"data$1\_size":$_ } @gsl_params;
	  for (@xs_args) {
	    if ($_->{name} =~ /^(n[1-9]?|w?stride[1-9]?)$/) {
	      #warn "XXXXXXXXXX '$_->{name}'\n";
	    }
	    else {
	      push @xs_args_censored, $_;
	    }
	  }
	}
	else {
	  @xs_args_censored = @xs_args;
	}
	
        my @if0 = map { "!$_->{name}" } @init_arr;
	my @if1 = map { "$_->{name}size<0" } @init_arr;
	my @if2 = map { "$_(\%d)" } @init_arr_size_if;
	my @if3 = map { "$n_first!=$_" } @n_other;
	push(@allfnc, { orig => $orig, svfunc => $svfunc, gsl_name => $gsl_name, xs_name => $xs_name,
			r2val => $r2val,
			result_array_size_msg => $result_array_size_msg,
			result_array_size => $result_array_size,
			xs_params => join(',',(map { $_->{name} } @xs_args_censored)), 
			xs_partypes => \@xs_partypes,
			xs_args => \@xs_args_censored,
			init_decl => \@init_decl,
			init_arr => \@init_arr,
			init_arr_size => \@init_arr_size,
			init_arr_size_if => join(' || ', @if1),
			init_arr_size_if_n => join(' || ', @if3),
			init_arr_alloc_if => join(' || ', @if0),
			init_arr_size_warn1 => join(', ', @if2),
			init_arr_size_warn2 => join(", ", @init_arr_size_if),
			gsl_params => join(',',@gsl_params),
		      } );
      }
      else {
        die "invalid line '$_'" unless m|^//|;
      }
    }

    my $tt = Template->new();
    #warn Dumper(\@allfnc);
    warn "tt.start\n";
    my $rv = $tt->process($template, { allfnc => \@allfnc, intro => $intro } , $outxs) or die "ERROR during process()\n";
    warn "tt.done\n";
  }  
  else {
    die "ERROR invalid file '$list'";
  }
}

proc_list('Math::EasyGSL', 'functions.list', qr/^#/, 'Functions.tt', rel2abs('..\lib\Math\EasyGSL\Functions.xs'));
proc_list('Math::EasyGSL::PDF', 'pdf.list', qr/^#/, 'PDF.tt', rel2abs('..\lib\Math\EasyGSL\PDF.xs'));
proc_list('Math::EasyGSL::CDF', 'cdf.list', qr/^#/, 'CDF.tt', rel2abs('..\lib\Math\EasyGSL\CDF.xs'));
proc_list('Math::EasyGSL::Random', 'random.list', qr/^#/, 'Random.tt', rel2abs('..\lib\Math\EasyGSL\Random.xs'));
proc_list('Math::EasyGSL::Statistics', 'statistics.list', qr/^#/, 'Statistics.tt', rel2abs('..\lib\Math\EasyGSL\Statistics.xs'));

#proc_list(undef, 'pdf.list', qr/^\?/, 'PDF.tt', 'PDF_questionmark.xs');
#proc_list(undef, 'random.list', qr/^\?/, 'Random.tt', 'Random_questionmark.xs');
#proc_list(undef, 'statistics.list', qr/^\?/, 'Statistics.tt', 'Statistics_questionmark.xs');
#proc_list(undef, 'statistics.list', qr/^\?/, 'Statistics.tt', 'Statistics_questionmark.xs');
