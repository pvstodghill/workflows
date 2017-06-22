#! /usr/bin/env perl

use strict;
use warnings FATAL => 'all';
use Carp::Always;

use FindBin;
use lib "$FindBin::Bin";
use Pvs::Psm;
use Pvs::Misc;

use constant { TRUE => 1, FALSE => 0 };

# ------------------------------------------------------------------------
# Process command line
# ------------------------------------------------------------------------

use File::Basename;
use Getopt::Std;

my $progname = basename($0);
my $getopt_usage = "Usage: cat psms.tsv | $progname [options] proteins.faa annotation.gff\n";
my $getopt_args = "";

our $opt_h; $getopt_args .= 'h';
$getopt_usage .= "-h - print help\n";

our $opt_p = 2; $getopt_args .= 'p:';
$getopt_usage .= "-p N - print proteins supported by >= N peptides [$opt_p]\n";

our $opt_o = FALSE; $getopt_args .= 'o';
$getopt_usage .= "-o - generate 'supported orfs'\n";

our $opt_r = FALSE; $getopt_args .= 'r';
$getopt_usage .= "-r - generate 'supported regions'\n";

our $opt_s = "unknown"; $getopt_args .= 's:';
$getopt_usage .= "-s STR - use STR for \"source\" field. [$opt_s]\n";

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

my ($proteins_faa,$annotation_gff,$ignored) = @ARGV;

if ( !defined($annotation_gff) || defined($ignored) ) {
  usage(1);
}

if ( !$opt_o && !$opt_r ) {
  print STDERR "Either -o or -r is required.\n";
  exit(1);
}

# ------------------------------------------------------------------------
# Read and store $proteins_faa
# ------------------------------------------------------------------------

if ( $opt_v ) {
  print STDERR "# $progname: Reading fasta file $proteins_faa\n";
}

my %metadata_of;

open(my $faa_fh,"<",$proteins_faa)
  || die "Cannot open <<$proteins_faa>>,";
while (<$faa_fh>) {
  chomp;
  if ($_ eq "") {
    next;
  }
  if ($_ !~ /^>/) {
    next;
  }
  my ($faa_protein,$faa_metadata) = parse_metadata($_);
  $metadata_of{$faa_protein} = $faa_metadata;
}
close($faa_fh);

# ------------------------------------------------------------------------
# Read and store the PSM's by protein id
# ------------------------------------------------------------------------

if ( $opt_v ) {
  print STDERR "# $progname: Reading features from $annotation_gff\n";
}

# store and retrieve half closed intervals [lb,ub)
use Set::IntervalTree;

my $spatial_indices = {};

open(my $fh, "<", $annotation_gff) || die "Cannot open <<$annotation_gff>>,";
while (<$fh>) {
  if ( /^#/ ) {
    next;
  }
  chomp;
  my ($seqname,$source,$feature,$start,$end,
      $score,$strand,$frame,$attributes) = split(/\t/,$_);

  my $entry = {
	       seqname => $seqname,
	       source => $source,
	       feature => $feature,
	       start => $start,
	       end => $end,
	       score => $score,
	       strand => $strand,
	       frame => $frame,
	       attributes => $attributes
	      };

  my $accession = $seqname;
  $accession =~ s/\.[0-9]+$//;

  ($start < $end) || die "<<$_>>,";
  my $spatial_index = $spatial_indices->{$accession};
  if (!$spatial_index) {
    $spatial_index = $spatial_indices->{$accession} =
      Set::IntervalTree->new;
  }
  $spatial_index->insert($entry,$start,$end+1);
}
close $fh;

# ------------------------------------------------------------------------
# Read and store the PSM's by protein id
# ------------------------------------------------------------------------

if ( $opt_v ) {
  print STDERR "# $progname: Reading PSM's from stdin\n";
}

my %psms_of;
while (my $psm = read_psm_v2()) {
  my $protein = $psm->{mapping}->{protein};
  if (!defined($psms_of{$protein})) {
    $psms_of{$protein} = [];
  }
  push @{$psms_of{$protein}}, $psm;
}

# ------------------------------------------------------------------------
# Filter and print protein features
# ------------------------------------------------------------------------

if ( $opt_v ) {
  print STDERR "# $progname: Printing GFF to stdout\n";
}

print "##gff-version 3\n";

my $source = $opt_s;
my $frame = ".";

# Deterministic output is nice.
foreach my $protein ( sort(keys(%psms_of)) ) {

  my $metadata = $metadata_of{$protein};
  my $accession = $metadata->{accession};
  my $version = $metadata->{version};
  my $orf_start = $metadata->{start};
  my $orf_end = $metadata->{end};
  my $strand = $metadata->{strand};

  my $protein_frame_index = frame_index($metadata);

  my $seqname = $accession;
  if (defined($version) && $version ne "") {
    $seqname .= ".".$version;
  }

  my @attributes;
  push @attributes, "orf_id=".$protein;

  my @psms = @{$psms_of{$protein}};
  (scalar(@psms) > 0) || die;

  if ( (scalar(@psms) < $opt_p) ) {
    next; 
  }
  push @attributes, "support=".scalar(@psms);

  my $score = 1.0;
  foreach my $psm ( @psms ) {
    $score *= $psm->{scores}->{pep}; #fixme: pep shouldn't be hardcoded
  }
  push @attributes, "pvalue=".$score;

  ($orf_start < $orf_end) || die;
  my $aref = $spatial_indices->{$accession}->fetch($orf_start,$orf_end+1);
  my @overlapping_annotations = @$aref;
  my @in_frame_annotations;
  foreach my $entry (@overlapping_annotations) {
    my $entry_frame_index = frame_index($entry);
    if ( $entry_frame_index != $protein_frame_index ) {
      next;
    }
    push @in_frame_annotations, $entry;
  }

  my $missing_cds = TRUE;
  # there can be multiple CDS's per ORF. E.g., Cmm382 CMM_2309 and CMM_PS_21
  my @cds_entries;
  my $pseudogene = FALSE;
  foreach my $entry (@in_frame_annotations) {
    if ( $entry->{feature} eq "CDS" ) {
      $missing_cds = FALSE;
      push @cds_entries, $entry;
    }
    if ( $entry->{attributes} =~ /;gene_biotype=pseudogene;/ ) {
      $pseudogene = TRUE;
    }
  }

  my $reannotate_cds = FALSE;
  if ( $missing_cds ) {
    push @attributes, "missing_cds";
  } else {
    foreach my $cds_entry ( @cds_entries ) {
      my $cds_start = $cds_entry->{start};
      my $cds_end = $cds_entry->{end};
      foreach my $psm ( @psms ) {
	my ($psm_start,$psm_end,$psm_strand) =
	  xlate_aa_coords($metadata,
			  $psm->{mapping}->{offset},
			  length($psm->{peptide}));
	if ( !($cds_start <= $psm_start && $psm_end <= $cds_end) ) {
	  $reannotate_cds = TRUE;
	}
      }
    }
  }

  if ( $pseudogene ) {
    push @attributes, "pseudogene";
  }
  if ( $reannotate_cds ) {
    push @attributes, "reannotate_cds";
  }

  my $attributes = join("; ",@attributes);

  if ( $opt_o ) {
    my $orf_feature = "ORF";
    print join("\t",$seqname,$source,$orf_feature,$orf_start,$orf_end,
	       $score,$strand,$frame,$attributes),"\n";
  }
  if ( $opt_r ) {
    my $orf_length = ($orf_end - $orf_start + 1) / 3;
    my $region_feature = "CDS_region";
    my $min_offset = min(map {$_->{mapping}->{offset}} @psms);
    my ($region_start,$region_end,$region_strand) =
      xlate_aa_coords($metadata,
		      $min_offset,
		      $orf_length - $min_offset);
    print join("\t",$seqname,$source,$region_feature,$region_start,$region_end,
	       $score,$strand,$frame,$attributes),"\n";
  }
}

sub frame_index {
  my ($x) = @_;
  my $strand = $x->{strand};
  my $start = $x->{start};
  my $end = $x->{end};
  if ($strand eq "+") {
    return  $start%3+1;
  } else {
    return -($end%3+1);
  }
}

# ------------------------------------------------------------------------
# Done
# ------------------------------------------------------------------------
