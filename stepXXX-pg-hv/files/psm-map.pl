#! /usr/bin/env perl

use strict;
use warnings FATAL => 'all';
use Carp::Always;

use FindBin;
use lib "$FindBin::Bin";
use Pvs::Psm;

use constant { TRUE => 1, FALSE => 0 };

# ------------------------------------------------------------------------
# Process command line
# ------------------------------------------------------------------------

use File::Basename;
use Getopt::Std;

my $progname = basename($0);
my $getopt_usage = "Usage: cat psms.tsv | $progname [options] peptides.faa\n";
my $getopt_args = "";

our $opt_S = 1;
$getopt_usage .= "-S N - random seed, for -p 2 [$opt_S]\n";
$getopt_args .= 'S:';

our $opt_h;
$getopt_usage .= "-h - print help\n";
$getopt_args .= 'h';

our $opt_n = 5;
$getopt_usage .= "-n N - hash string size [$opt_n]\n";
$getopt_args .= 'n:';

our $opt_p = 1;
$getopt_usage .= "-p 1 - Only report peptides with unique genomic hits (default)\n";
$getopt_usage .= "-p 2 - Report one random genomic hit per peptide\n";
$getopt_usage .= "-p 3 - Report all genomic hits\n";
$getopt_args .= 'p:';

our $opt_v;
$getopt_usage .= "-v - run verbosely\n";
$getopt_args .= 'v';

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

my ($faa_file,$ignored) = @ARGV;

if ( !defined($faa_file) || defined($ignored) ) {
  usage(1);
}

if ( $opt_p ne "1" && $opt_p ne "2" && $opt_p ne "3" ) {
  usage(1);
}

# ------------------------------------------------------------------------
# shuffle lists
# ------------------------------------------------------------------------

srand($opt_S);

sub shuffle {
  my @l = @_;
  @l = map { [rand(),$_] } @l;
  @l = sort { $a->[0] cmp $b->[0] } @l;
  @l = map { $_->[1] } @l;
  return @l;
}

# ------------------------------------------------------------------------
# Read and store $faa_file
# ------------------------------------------------------------------------

if ( $opt_v ) {
  print STDERR "# $progname: Reading fasta file\n";
}

# my %metadata_of;

my $faa_defline;
my $faa_sequence;

init_buckets($opt_n);

sub store_protein {
  if ( !defined($faa_defline) ) {
    return;
  }

  my ($faa_protein,$faa_metadata) = parse_metadata($faa_defline);
  # $metadata_of{$faa_protein} = $faa_metadata;
  add_sequence($faa_protein, $faa_sequence);

  $faa_defline = $faa_sequence = undef;
}

open(my $faa_fh,"<",$faa_file)
  || die "Cannot open <<$faa_file>>,";
while (<$faa_fh>) {
  chomp;
  if ($_ eq "") {
    next;
  }
  if ($_ !~ /^>/) {
    $faa_sequence .= $_;
    next;
  }
  store_protein();
  $faa_defline = $_;
  $faa_sequence = "";
}
store_protein();
close($faa_fh);


# ------------------------------------------------------------------------
# Print bucket stats
# ------------------------------------------------------------------------

if ( $opt_v ) {
  Pvs::Psm::dump_bucket_stats();
}


# ------------------------------------------------------------------------
# Read and process $pep_xml_file
# ------------------------------------------------------------------------

if ( $opt_v ) {
  print STDERR "# $progname: Reading PSM's from stdin, printing mapped PSM's to stdout\n";
}

while (my $psm = read_psm_v2()) {
  map_and_print_psm($psm);
}

# ------------------------------------------------------------------------
# map psm's to the proteome and print new psms
# ------------------------------------------------------------------------

sub map_and_print_psm {
  my ($psm) = @_;

  my $spectrum = $psm->{spectrum};
  my $peptide = $psm->{peptide};

  my @results = search_sequences($peptide);
  (scalar(@results) >= 1) || die;

  # sanity checks
  my $match_protein = $psm->{mapping}->{protein};
  if ( defined($match_protein) ) {
    my $found = FALSE;
    foreach my $e ( @results ) {
      my ($hit_protein,$offset) = @$e;
      if ( $hit_protein eq $match_protein ) {
	$found = TRUE;
      }
    }
    ( $found ) || die "Not found - <<$spectrum>><<$peptide>>,";
  }
  
  my $num_locations = scalar(@results);

  my @reported;
  if ( $opt_p == 1 ) {
    if ( $num_locations > 1 ) {
      return;
    }
    @reported = @results;
  } elsif ( $opt_p == 2 ) {
    my @l = shuffle(@results);
    @reported = ($l[0]);
  } else {
    ($opt_p == 3) || die;
    @reported = @results;
  }

  my $psm2 = { %$psm };
  foreach my $e ( @reported ) {
    my ($hit_protein,$offset) = @$e;
    $psm2->{mapping} = {
			protein => $hit_protein,
			offset => $offset,
			length => length($peptide),
			num_locations => $num_locations
		      };
    print_psm_v2($psm2);
  }

}

# ------------------------------------------------------------------------

if ( $opt_v ) {
  print STDERR "# $progname: Done.\n";
}
