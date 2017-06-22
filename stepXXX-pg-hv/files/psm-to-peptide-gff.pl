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

our $opt_a; $getopt_args .= 'a:';
$getopt_usage .= "-a STR - extra attributes\n";

our $opt_h; $getopt_args .= 'h';
$getopt_usage .= "-h - print help\n";

our $opt_m; $getopt_args .= 'm';
$getopt_usage .= "-m - just modifications, no peptides\n";

our $opt_p; $getopt_args .= 'p';
$getopt_usage .= "-p - just peptides, no modifications\n";

our $opt_s = "unknown"; $getopt_args .= 's:';
$getopt_usage .= "-s STR - use STR for \"source\" field. [$opt_s]\n";

our $opt_v; $getopt_args .= 'v';
$getopt_usage .= "-v - run verbosely\n";

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

if ( !defined($opt_m) &&  !defined($opt_p) ) {
  $opt_m = $opt_p = TRUE;
}

# ------------------------------------------------------------------------
# Read and store $faa_file
# ------------------------------------------------------------------------

if ( $opt_v ) {
  print STDERR "# $progname: Reading $faa_file\n";
}

my %metadata_of;
my %sequence_of;

my $faa_defline;
my $faa_sequence;

sub store_protein {
  if ( !defined($faa_defline) ) {
    return;
  }

  my ($faa_protein,$faa_metadata) = parse_metadata($faa_defline);
  $metadata_of{$faa_protein} = $faa_metadata;
  $sequence_of{$faa_protein} = $faa_sequence;

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
# Read and process the PSM's
# ------------------------------------------------------------------------

my $next_modification_colour = 2; # Start with red
my %modification_colours;

if ( $opt_v ) {
  print STDERR "# $progname: Reading PSM's from stdin, printing GFF to stdout.\n";
}

print "##gff-version 3\n";
while (my $psm = read_psm_v2()) {
  emit_gff_hit($psm);
}

# ------------------------------------------------------------------------

sub emit_gff_hit {
  my ($psm) = @_;

  my $sample = $psm->{sample};
  my $spectrum = $psm->{spectrum};
  my $peptide = $psm->{peptide};
  my $p_value = $psm->{scores}->{p_value};
  my $q_value = $psm->{scores}->{q_value};
  my $pep = $psm->{scores}->{pep};
  my $xcorr = $psm->{scores}->{xcorr};
  my $hit_protein = $psm->{mapping}->{protein};
  my $hit_offset = $psm->{mapping}->{offset};
  my $num_locations = $psm->{mapping}->{num_locations};

  my $protein_metadata = $metadata_of{$hit_protein};
  ( defined($protein_metadata) ) || die "Not found in DB: $hit_protein <<$peptide>>,";

  my $N_terminal_M = $protein_metadata->{N_terminal_M};
  my $N_terminal_MAP_processed = $protein_metadata->{N_terminal_MAP_processed};
  my $is_decoy = $protein_metadata->{is_decoy};
  my $accession = $protein_metadata->{accession};
  my $version = $protein_metadata->{version};

  my $seqname = $accession.".".$version;
  my $source = $opt_s;
  my $score = $xcorr;
  my $frame = ".";

  # ------------------------------------------------------------------------

  my @attributes_common = ();
  push @attributes_common, [ "sample", multiple_to_string($sample) ];
  push @attributes_common, [ "spectrum", multiple_to_string($spectrum) ];
  push @attributes_common,[ "peptide", $peptide ];
  push @attributes_common,[ "pvalue", $p_value ];
  push @attributes_common,[ "qvalue", $q_value ];
  push @attributes_common,[ "PEP", $pep ];
  push @attributes_common,[ "XCorr", $xcorr ];
  if ($is_decoy) {
    push @attributes_common, [ "decoy" ];
  }
  if ($hit_offset == 0) {
    push @attributes_common, [ "offsetzero" ];
  }
  if ($N_terminal_M) {
    push @attributes_common, [ "N_terminal_M" ] ;
  }
  if ($N_terminal_MAP_processed) {
    push @attributes_common, [ "N_terminal_MAP_processed" ] ;
  }

  # ------------------------------------------------------------------------

  if ( $opt_p ) {

    my $feature_pep = "misc_peptide";
    my ($start_pep, $end_pep, $strand_pep) =
      xlate_aa_coords($metadata_of{$hit_protein},
		      $hit_offset,
		      length($peptide));

    my @attributes_pep = ();

    push @attributes_pep, [ "locations", $num_locations ];
    my $colour = ( $num_locations == 1 ) ? 3 : 7;
    push @attributes_pep, [ "colour", $colour ];

    my $attributes_pep = unparse_attrs(@attributes_common, @attributes_pep);

    print join("\t",$seqname,$source,$feature_pep,$start_pep,$end_pep,
	       $score,$strand_pep,$frame,$attributes_pep),"\n";

  }

  # --------------------------------------------------

  if ( $opt_m ) {

    my $feature_mod = "modified_amino_acid_feature";

    foreach my $mod ( @{$psm->{modifications}} ) {
      my ($mod_offset,$mod_name) = @$mod;
      my ($start_mod, $end_mod, $strand_mod) =
	xlate_aa_coords($metadata_of{$hit_protein},
			$hit_offset+$mod_offset,
			1);
      my @attributes_mod = ();
      push @attributes_mod,[ "modification", "\"".$mod_name."\"" ];
      my $colour = $modification_colours{$mod_name};
      if ( !defined($modification_colours{$mod_name}) ) {
	$colour = $modification_colours{$mod_name} = $next_modification_colour++;
      }
      push @attributes_mod, [ "colour", $colour ];

      my $attributes_mod = unparse_attrs(@attributes_common,
					 @attributes_mod);
      print join("\t",$seqname,$source,$feature_mod,$start_mod,$end_mod,
		 $score,$strand_mod,$frame,$attributes_mod),"\n";
    }

  }

}

sub multiple_to_string {
  my ($x) = @_;
  if ( !ref($x) ) {
    return $x;
  } else {
    return join(",",@$x);
  }
}

sub unparse_attrs {
  return join(";", (map { unparse_attr(@{$_}) } @_));
}

sub unparse_attr {
  my ($name,$val) = @_;
  (defined($name)) || die;
  my $attr = $name;
  if ( defined($val) ) {
    $attr .= '="'.$val.'"';
  }
  return $attr;
}

# ------------------------------------------------------------------------

if ( $opt_v ) {
  print STDERR "# $progname: Done.\n";
}
