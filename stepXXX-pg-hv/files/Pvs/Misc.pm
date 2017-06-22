package Pvs::Misc;

use Exporter qw(import);

our @ISA    = qw(Exporter);
our @EXPORT = qw(
		  unique
		  min
		  max
		  sum
	       );
our @EXPORT_OK = qw(
		  );

use strict;
use warnings FATAL => 'all';
use Carp::Always;

use constant { TRUE => 1, FALSE => 0 };

sub unique {
  my ($cmp,@in) = @_;
  if ( scalar(@in) == 0 ) {
    return ();
  }
  my ($x,@sorted) = sort { &$cmp($a,$b) } @in;
  my @out = ($x);
  foreach my $y ( @sorted ) {
    my $x = $out[-1];
    $a = $x; $b = $y;
    if ( &$cmp($x,$y) != 0 ) {
      push @out, $y;
    }
  }
  return @out;
}

sub min {
  my ($x, @l) = @_;
  foreach my $y ( @l ) {
    if ( $x > $y ) {
      $x = $y;
    }
  }
  return $x;
}

sub max {
  my ($x, @l) = @_;
  foreach my $y ( @l ) {
    if ( $x < $y ) {
      $x = $y;
    }
  }
  return $x;
}

sub sum {
  my $sum = 0;
  foreach my $x ( @_ ) {
    $sum += $x;
    return $sum;
  }
}

1;
