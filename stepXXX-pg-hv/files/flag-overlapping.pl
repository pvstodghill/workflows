#! /usr/bin/env perl

use strict;
use warnings FATAL => 'all';
use Carp::Always;

# use FindBin;
# use lib "$FindBin::Bin";
# use Xyzzy;

use constant { TRUE => 1, FALSE => 0 };

# ------------------------------------------------------------------------
# Process command line
# ------------------------------------------------------------------------

use File::Basename;
use Getopt::Std;

my $progname = basename($0);
my $getopt_usage = "Usage: cat features.gff | $progname [options]\n";
my $getopt_args = "";

our $opt_a = "overlaps"; $getopt_args .= 'a:';
$getopt_usage .= "-a STR - attribute to add for overlapps [\"$opt_a\"]\n";

our $opt_h; $getopt_args .= 'h';
$getopt_usage .= "-h - print help\n";

our $opt_o; $getopt_args .= 'o';
$getopt_usage .= "-o - only print overlapping features\n";

our $opt_v; $getopt_args .= 'v';
$getopt_usage .= "-v - be verbose\n";

sub usage {
  print STDERR $getopt_usage;
  exit(@_);
}

my $stat = getopts($getopt_args);
if (!$stat) {
  usage(1);
}
if ($opt_h) {
  usage();
}

# ------------------------------------------------------------------------
# Read input GFF
# ------------------------------------------------------------------------

if ( $opt_v ) {
  print STDERR "# $progname: Reading GFF features from stdin\n";
}

# store and retrieve half closed intervals [lb,ub)
use Set::IntervalTree;

my $spatial_indices = {};

my @entries;

my $n = 0;
while (<STDIN>) {
  if ( /^#/ ) { next; }
  chomp;
  my ($seqname,$source,$feature,$start,$end,
    $score,$strand,$frame,$attributes) = split(/\t/,$_);

  my $accession = $seqname;
  $accession =~ s/\.[0-9]+$//;

  my $entry = {
	       id => $n++,
	       seqname => $seqname,
	       accession => $accession,
	       source => $source,
	       feature => $feature,
	       start => $start,
	       end => $end,
	       score => $score,
	       strand => $strand,
	       frame => $frame,
	       attributes => $attributes,
	       overlap => FALSE,
	       overlap_opposite => FALSE
	       };

  push @entries, $entry;

  ($start < $end) || die "<<$_>>,";
  my $spatial_index = $spatial_indices->{$accession};
  if (!$spatial_index) {
    $spatial_index = $spatial_indices->{$accession} =
      Set::IntervalTree->new;
  }
  $spatial_index->insert($entry,$start,$end+1);
}

# ------------------------------------------------------------------------
# Compute overlaps
# ------------------------------------------------------------------------

if ( $opt_v ) {
  print STDERR "# $progname: Computing overlapping features\n";
}

foreach my $entry ( @entries ) {
  my $accession = $entry->{accession};
  my $start = $entry->{start};
  my $end = $entry->{end};
  my $aref = $spatial_indices->{$accession}->fetch($start,$end+1);
  my @overlapping_entries = @$aref;
  # an entry will overlap with itself
  if ( scalar(@overlapping_entries) == 1 ) { next; }
  foreach my $entry2 ( @overlapping_entries ) {
    $entry2->{overlap} = TRUE;
    if ( $entry->{strand} ne $entry2->{strand} ) {
      $entry->{overlap_opposite} = $entry2->{overlap_opposite} = TRUE;
    }
  }
}


# ------------------------------------------------------------------------
# Print results
# ------------------------------------------------------------------------

if ( $opt_v ) {
  print STDERR "# $progname: Printing output GFF \n";
}

print "##gff-version 3\n";

foreach my $entry ( @entries ) {

  my $seqname = $entry->{seqname};
  my $source = $entry->{source};
  my $feature = $entry->{feature};
  my $start = $entry->{start};
  my $end = $entry->{end};
  my $score = $entry->{score};
  my $strand = $entry->{strand};
  my $frame = $entry->{frame};
  my $attributes = $entry->{attributes};

  my $overlap = $entry->{overlap};
  my $overlap_opposite = $entry->{overlap_opposite};

  if ($opt_o && !$overlap) {
    next;
  }

  if ( $overlap ) {
    if ( $attributes eq "" ) {
      $attributes = $opt_a;
    } else {
      $attributes .= ";".$opt_a;
    }
  }

  if ( $overlap_opposite ) {
    if ( $attributes eq "" ) {
      $attributes = $opt_a."_opposite";
    } else {
      $attributes .= ";".$opt_a."_opposite";
    }
  }

  print join("\t",$seqname,$source,$feature,$start,$end,
	     $score,$strand,$frame,$attributes),"\n";
}
