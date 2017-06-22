#! /usr/bin/env perl

use strict;
use warnings FATAL => 'all';
use Carp::Always;

use FindBin;
use lib "$FindBin::Bin";
use Pvs::Psm;

use constant { TRUE => 1, FALSE => 0 };

# ------------------------------------------------------------------------
# Process the command line
# ------------------------------------------------------------------------

use File::Basename;
use Getopt::Std;

my $progname = basename($0);

our $opt_c; # allowable number of errors in final set < $opt_c
our $opt_h;
our $opt_i; # attribute providing p-value
our $opt_m; # p.adjust method
our $opt_o; # attribute receiving q-value
our $opt_v;

sub usage {
  print STDERR "Usage: $progname [options] ...\n";
  print STDERR "-c CUTOFF - allowable number of errors in final set < CUTOFF \n";
  print STDERR "-h - print help\n";
  print STDERR "-i NAME - attribute to use as input q-value\n";
  print STDERR "-v - be verbose\n";
  print STDERR "Required: -c -i. \n";
  exit(@_);
}

my $stat = getopts('c:hi:v');
if (!$stat) {
  usage(1);
}
if ($opt_h) {
  usage();
}

if ( !defined($opt_c) || !defined($opt_i) ) {
  usage(1);
}


# ------------------------------------------------------------------------
# Read, apply filters, print
# ------------------------------------------------------------------------

if ( $opt_v ) {
  print STDERR "# $progname: Reading PSM's from stdin, printing filtered PSM's to stdout\n";
}

my @psms;
while (my $psm = read_psm_v2()) {
  push @psms, $psm;
  my $qvalue = $psm->{scores}->{$opt_i};
  if ( !defined($qvalue) ) {
    print STDERR "No value for attribute <<<$opt_i>>>, entry:\n";
    print STDERR Dumper($psm),"\n";
    exit(1);
  }  
}

# ------------------------------------------------------------------------
# Bin the PSM's by FDR value
# ------------------------------------------------------------------------

my %bins;
foreach my $psm ( @psms ) {
  my $fdr = $psm->{scores}->{$opt_i};
  my $a = $bins{$fdr};
  if (!defined($a)) {
    $a = $bins{$fdr} = [];
  }
  push @$a, $psm;
}

# ------------------------------------------------------------------------
# Output PSM bin's until we meet the $opt_c cutoff.
# ------------------------------------------------------------------------

my $count = 0;
foreach my $fdr ( sort {$a<=>$b} (keys %bins) ) {
  my $r = $bins{$fdr};
  $count += scalar(@$r);
  if ( $fdr * $count >= $opt_c ) {
    last;
  }
  # deterministic resuls are nice
  foreach my $psm (sort {cmp_psms($a,$b)} (@$r)) {
    print_psm_v2($psm);
  }
}
