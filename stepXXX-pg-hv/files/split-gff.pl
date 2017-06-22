#! /usr/bin/env perl

use strict;
use warnings;

use constant { TRUE => 1, FALSE => 0 };

# Script to demux a GFF file into accession specific GFF files.

# Usage:
# $ rm -f *.gff
# $ cat ..../input.gff | split-gff
# $ ls
# AE016853.gff  CP000058.gff  CP000075.gff
# $ rm -f *.gff
# $ cat ..../input.gff | split-gff foo
# $ ls
# AE016853_foo.gff  CP000058_foo.gff  CP000075_foo.gff
# $ 

# ------------------------------------------------------------------------

use File::Basename;
use Getopt::Std;

our $opt_d = ".";
our $opt_h;
our $opt_n;

sub usage {
  my $progname = basename($0);
  print STDERR "Usage: $progname [options] [tag]\n";
  print STDERR "-d DIR - output to DIR [.]\n";
  print STDERR "-h - print help\n";
  print STDERR "-n - strip version from accession\n";
  exit(@_);
}

my $stat = getopts('d:hn');
if (!$stat) {
  usage(1);
}
if ($opt_h) {
  usage();
}

my ($tag) = @ARGV;
if (defined($tag)) {
  $tag = "_".$tag;
} else {
  $tag = "";
}

# ------------------------------------------------------------------------

my $gff3 = FALSE;

my %seqnames;

while (<STDIN>) {
  chomp;
  $_ =~ s/\r+//;
  if ( /^#/ ) {
    if ( $_ eq "##gff-version 3" ) {
      $gff3 = TRUE;
    }
    if ( $_ eq "##FASTA" ) {
      last;
    }
    next;
  }
  my ($seqname) = split(/\t/,$_);
  if ( $opt_n ) {
    $seqname =~ s/\.[0-9]+$//;
  }
  ( $seqname !~ /\|/ ) || die "cannot parse seqname=<<$seqname>>,";
  if (!defined($seqname)) {
    $seqnames{$seqname} = [];
  }
  push @{$seqnames{$seqname}}, $_;
}

# ------------------------------------------------------------------------

foreach my $seqname ( keys %seqnames ) {
  my $out_name =  $opt_d."/".$seqname.$tag.".gff";
  open(F, ">", $out_name) || die "Cannot open for writing <<$out_name>>,";
  if ( $gff3 ) {
    print F "##gff-version 3\n";
  }
  foreach my $l ( @{$seqnames{$seqname}} ) {
    print F $l,"\n";
  }
  close F;
}
