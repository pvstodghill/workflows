#! /usr/bin/env perl

use strict;
use warnings;
use Carp::Always;

use constant { TRUE => 1, FALSE => 0 };

# ------------------------------------------------------------------------
# Process command line
# ------------------------------------------------------------------------

use File::Basename;
use Getopt::Std;

our $opt_h;

sub usage {
  my $progname = basename($0);
  print STDERR "Usage: $progname [options] locus_names.gff regulon.gff > regulon.txt \n";
  print STDERR "-h - print help\n";
  exit(@_);
}

my $stat = getopts('h');
if (!$stat) {
  usage(1);
}
if ($opt_h) {
  usage();
}

if ( scalar(@ARGV) != 2 ) {
  usage(1);
}

my ($locus_gff,$regulon_gff) = @ARGV;

# ------------------------------------------------------------------------
# Read the locus gff file
# ------------------------------------------------------------------------

my %name_of_locus;

open(F, "<", $locus_gff) || die;
while (<F>) {
  if ( /^#/) { next; }
  chomp;
  my ($seqname,$source,$feature,$start,$end,
    $score,$strand,$frame,$attributes) = split(/\t/,$_);
  my $locus = "$start..$end/$strand";
  if ($feature eq "gene") {
    ( $attributes =~ /(^|;)locus_tag=([^;]+)/ ) || die;
    my $name = $2;
    $name_of_locus{$locus} = $name;
  } else {
    die "Non-gene feature in locus file <<$_>>,";
  }
}
close F;


# ------------------------------------------------------------------------
# Read the regulon table;
# ------------------------------------------------------------------------

my %printed_name;

open(F, "<", $regulon_gff) || die;
while (<F>) {
  if ( /^#/) { next; }
  chomp;
  my ($seqname,$source,$feature,$start,$end,
    $score,$strand,$frame,$attributes) = split(/\t/,$_);
  my $locus = "$start..$end/$strand";
  my $name = $name_of_locus{$locus};
  if ( !defined($name) ) {
    die "Cannot find locus name, <<$_>>,";
  } elsif ( $printed_name{$name} ) {
    die "Printed locus name twice?, <<$locus>> <<$name>>,";
  }
  print "$name\n";
}
close F;
