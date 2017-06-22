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
my $getopt_usage = "Usage: $progname [options] file1.pep.xml file2.pep.xml ...\n";
my $getopt_args = "";

our $opt_h; $getopt_args .= 'h';
$getopt_usage .= "-h - print help\n";

our $opt_l = "unspecified"; $getopt_args .= 'l:';
$getopt_usage .= "-l LABEL - sample label [$opt_l]\n";

our $opt_m = 1; $getopt_args .= 'm:';
$getopt_usage .= "-m 1 - Report best <search_hit> for each <spectrum_query> (default)\n";
$getopt_usage .= "-m 2 - Report all <search_hit>'s for each <spectrum_query>\n";

our $opt_q = 1e-3; $getopt_args .= 'q:';
$getopt_usage .= "-q NUM - q-value cutoff [$opt_q]\n";

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

if ( scalar(@ARGV) == 0 ) {
  usage(1);
}

if ( $opt_m ne "1" && $opt_m ne "2" ) {
  usage(1);
}

# ------------------------------------------------------------------------
# Parser/Printer routines
# ------------------------------------------------------------------------

use XML::LibXML::Reader;

sub process_pep_xml_file {
  my ($pep_xml_file) = @_;
  my $reader = XML::LibXML::Reader->new(location => $pep_xml_file);
  ($reader) || die;
  reset_parse_state();
  while ($reader->read()) {
    if ($reader->nodeType == XML_READER_TYPE_ELEMENT) {
      if (FALSE) {
      } elsif ($reader->localName() eq "parameter") {
	bo_parameter(get_attrs($reader));
      } elsif ($reader->localName() eq "spectrum_query") {
	bo_spectrum_query(get_attrs($reader));
      } elsif ($reader->localName() eq "search_hit") {
	bo_search_hit(get_attrs($reader));
      } elsif ($reader->localName() eq "search_score") {
	bo_search_score(get_attrs($reader));
      } elsif ($reader->localName() eq "percolator_result") {
	bo_percolator_result(get_attrs($reader));
      } elsif ($reader->localName() eq "mod_aminoacid_mass") {
	bo_mod_aminoacid_mass(get_attrs($reader),$reader->readOuterXml());
      }
    } elsif ($reader->nodeType == XML_READER_TYPE_END_ELEMENT) {
      if (FALSE) {
      } elsif ($reader->localName() eq "percolator_result") {
	;
      } elsif ($reader->localName() eq "search_hit") {
	eo_search_hit(get_attrs($reader));
      } elsif ($reader->localName() eq "spectrum_query") {
	eo_spectrum_query(get_attrs($reader));
      }
    }
  }
  $reader->close();
}

sub get_attrs {
  my ($reader) = @_;
  my $h = {};
  while ($reader->moveToNextAttribute()) {
    my $attr = $reader->localName();
    my $value = $reader->value();
    $h->{$attr} = $value;
  }
  return $h;
}

# ------------------------------------------------------------------------

my @mods;

my $spectrum_query;
my @spectrum_hits;
my $search_hit;
my $percolator_result;
my $xcorr_score;
my @modifications;

sub reset_parse_state {
  @mods = ();
 $spectrum_query = undef;
 @spectrum_hits = ();
 $search_hit = undef;
 $percolator_result = undef;
 $xcorr_score = undef;
 @modifications = ();
}

# ------------------------------------------------------------------------

sub bo_parameter {
  my ($attrs) = @_;
  my $v = $attrs->{value};
  if (!$v) {
  } elsif ( $v =~ /^(.*) \/ (.*) Da \((.*)\)$/) {
    # PD 1.4 "^name / delta Da (aas)$"
    push @mods, [$1, $3, $2];
  } elsif ( $v =~ /^<Modification .+ AminoAcids="([^\"]*)" Name="([^\"]*)" .+ DeltaAverageMass="([^\"]*)" .* \/>$/ ) {
    # PD 2 "<Modification .* AminoAcids="aas" Name="name" .* DeltaMass="delta" .* />"
    push @mods, [$2, $1, $3];
  }

}

sub bo_spectrum_query {
  my ($attrs) = @_;
  (!defined($spectrum_query)) || die;
  $spectrum_query = $attrs;
}

sub eo_spectrum_query {
  my ($attrs) = @_;
  emit_spectrum_hits($spectrum_query,@spectrum_hits);
  $spectrum_query = undef;
  @spectrum_hits = ();
}
sub bo_search_hit {
  my ($attrs) = @_;
  (!defined($search_hit)) || die ("spectrum=".$spectrum_query->{spectrum}.",");
  $search_hit = $attrs;
}

sub eo_search_hit {
  my ($attrs) = @_;
  if ( $percolator_result->{"q-Value"} <= $opt_q ) {
    push @spectrum_hits, [$search_hit, $percolator_result, $xcorr_score, @modifications];
  }
  $xcorr_score = undef;
  $percolator_result = undef;
  $search_hit = undef;
  @modifications = ();
}

sub bo_search_score {
  my ($attrs) = @_;
  if ( defined($attrs->{name}) && $attrs->{name} eq "XCorr" ) {
    (!defined($xcorr_score)) || die;
    $xcorr_score = $attrs->{value};
  }
}

sub bo_percolator_result {
  my ($attrs) = @_;
  (!defined($percolator_result)) || die;
  $percolator_result = $attrs;
}

my $unmodified_mass =
  {
   G => 57.021464,
   A => 71.037114,
   S => 87.032028,
   P => 97.052764,
   V => 99.068414,
   T => 101.047678,
   C => 103.009184,
   I => 113.084064,
   L => 113.084064,
   N => 114.042927,
   D => 115.026943,
   Q => 128.058578,
   K => 128.094963,
   E => 129.042593,
   M => 131.040485,
   H => 137.058912,
   F => 147.068414,
   R => 156.101111,
   Y => 163.063329,
   W => 186.079313
  };

my $mod_tolerance = 0.01;


sub bo_mod_aminoacid_mass {
  my ($attrs,$xml) = @_;

  my $peptide = $search_hit->{peptide};

  my $pos = $attrs->{position};
  my $mod_mass = $attrs->{mass};
  my $aa = substr($peptide,$pos-1,1);
  my $unmod_mass = $unmodified_mass->{$aa};
  my $observed_delta = abs($mod_mass - $unmod_mass);
  my @candidate_mods;
  foreach my $mod ( @mods ) {
    my ($name,$aas,$predicted_delta) = @$mod;
    my $x = abs($observed_delta - $predicted_delta);
    if ( $x <= $mod_tolerance ) {
      push @candidate_mods, [$pos-1,$name];
    }
  }
  if (scalar(@candidate_mods) != 1) {
    print "pos = $pos\n";
    print $xml,"\n";
    foreach my $mod ( @mods ) {
      my ($pos,$name) = @$mod;
      print "    $name\n";
    }
    die;
  }
  push @modifications, @candidate_mods;
}

# ------------------------------------------------------------------------

sub emit_spectrum_hits {
  my ($spectrum_query,@spectrum_hits) = @_;
  if ( scalar(@spectrum_hits) == 0 ) {
    return;
  }
  if ($opt_m == 1) {
    # descreasing sort on xcorr
    @spectrum_hits = sort { $b->[2] <=> $a->[2] } @spectrum_hits;
    my $spectrum_hit = $spectrum_hits[0];
    emit_spectrum_hit($spectrum_query, $spectrum_hit);
  } else {
    foreach my $spectrum_hit (@spectrum_hits) {
      emit_spectrum_hit($spectrum_query, $spectrum_hit);
    }
  }
}

sub emit_spectrum_hit {
  my ($spectrum_query,$spectrum_hit) = @_;
  my ($search_hit, $percolator_result, $xcorr_score,@modifications) = @$spectrum_hit;
  print_psm_v2({
		sample => $opt_l,
		spectrum => $spectrum_query->{spectrum},
		peptide => $search_hit->{peptide},
		modifications => [@modifications],
		scores => {
			   p_value => $percolator_result->{"probability"},
			   q_value => $percolator_result->{"q-Value"},
			   pep => $percolator_result->{"PEP"},
			   xcorr => $xcorr_score,
			  },
		mapping => {
			    protein => $search_hit->{protein},
			   }
	       });
}


# ------------------------------------------------------------------------
# Read and process $pep_xml_file
# ------------------------------------------------------------------------

foreach my $pep_xml_file ( @ARGV ) {
  if ( $opt_v ) {
    print STDERR "# $progname: Processing $pep_xml_file\n";
  }
  process_pep_xml_file($pep_xml_file);
}

# ------------------------------------------------------------------------
# Done.
# ------------------------------------------------------------------------

if ( $opt_v ) {
  print STDERR "# $progname: Done.\n";
}
