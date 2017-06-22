package Pvs::Psm;

use Exporter qw(import);

our @ISA    = qw(Exporter);
our @EXPORT = qw(

		  dump_bucket_stats
		  unparse_modifications
		  add_sequence
		  cmp_psms
		  init_buckets
		  parse_metadata
		  print_psm_v2
		  read_psm_v2
		  search_sequences
		  xlate_aa_coords
);
our @EXPORT_OK = qw(

);

use strict;
use warnings FATAL => 'all';
use Carp::Always;

use constant { TRUE => 1, FALSE => 0 };


# ------------------------------------------------------------------------
# Parse Dave's FAA meta-data
# ------------------------------------------------------------------------

sub parse_metadata {
  my ($defline) = @_;

  ($defline =~ /^>([^ ]+) +(.*)/) || die "defline=<<$defline>>,";
  my ($id,$raw) = ($1,$2);

  my $parsed = {};
  my @l1 = split(" ",$raw);

  $parsed->{N_terminal_M} = FALSE;
  $parsed->{N_terminal_MAP_processed} = FALSE;
  if ( $l1[0] eq "N-terminal_M" ) {
    $parsed->{N_terminal_M} = TRUE;
    shift @l1;
  } elsif ( $l1[0] eq "N-terminal_MAP_processed" ) {
    $parsed->{N_terminal_MAP_processed} = TRUE;
    shift @l1;
  }

  $parsed->{is_decoy} = FALSE;
  if ( $l1[0] eq "decoy" ) {
    $parsed->{is_decoy} = TRUE;
    shift @l1;
    shift @l1;
  }
  (scalar(@l1) == 2) || die "protein_metadata=<<$raw>>";

  my ($orf_id,$position) = @l1;
  my @l2 = split(/\./,$orf_id);
  (scalar(@l2) == 3) || die "orf_id=<<$orf_id>>,";
  my ($accession,$version,$orf_num) = @l2;
  my ($orf_start,$orf_end,$orf_strand);
  if ( $position =~ /^complement\(([0-9]+)\.\.([0-9]+)\)$/ ) {
    ($orf_start,$orf_end,$orf_strand) = ($1,$2,"-");
  } elsif ( $position =~ /^([0-9]+)\.\.([0-9]+)$/ ) {
    ($orf_start,$orf_end,$orf_strand) = ($1,$2,"+");
  } else {
    die "position = <<$position>>,";
  }
  ($orf_start <= $orf_end) || die "$orf_start > $orf_end,";

  $parsed->{id} = $id;
  $parsed->{accession} = $accession;
  $parsed->{version} = $version;
  $parsed->{start} = $orf_start;
  $parsed->{end} = $orf_end;
  $parsed->{strand} = $orf_strand;

  return ($id,$parsed);
}


# ------------------------------------------------------------------------
# Parsing and unparsing PSM files
# ------------------------------------------------------------------------

my $printed_header = FALSE;

sub print_psm_v2 {
  my ($psm) = @_;
  if (!$printed_header) {
    print join("\t",
	       "spectrum",
	       "peptide",
	       "sample",
	       "modifications",
	       "scores",
	       "mapping",
	      ),"\n";
    $printed_header = TRUE;
  }

  print join("\t",
	     unparse_multiple($psm->{spectrum}),
	     unparse_single($psm->{peptide}),
	     unparse_multiple($psm->{sample}),
	     unparse_modifications($psm->{modifications}),
	     unparse_named($psm->{scores}),
	     unparse_named($psm->{mapping}),
	    ),"\n";

};

sub unparse_modifications { return unparse_indexed(@_); }

# ------------------------------------------------------------------------

sub read_psm_v2 {
  my $l = <STDIN>;
  if (!$l) {
    return undef;
  }
  chomp $l;
  my ($spectrum,$peptide,$sample,$modifications,$scores,$mapping)
    = split("\t",$l);
  if ( $spectrum eq "spectrum" ) {
    return read_psm_v2();
  }
  return {
	  spectrum => parse_multiple($spectrum),
	  peptide => parse_single($peptide),
	  sample => parse_multiple($sample),
	  modifications => parse_indexed($modifications),
	  scores => parse_named($scores),
	  mapping => parse_named($mapping),
	 };

}

# ------------------------------------------------------------------------

sub parse_single {
  my ($s) = @_;
  my @entries = split(/; */,$s);
  ( scalar(@entries) == 1 ) || die "parse_single(<<$s>>),";
  return $s;
}


sub unparse_single {
  my ($x) = @_;
  (!ref($x)) || die "unparse_single(REF),";
  return $x;
}

# ------------------------------------------------------------------------

sub parse_multiple {
  my ($s) = @_;
  my @entries = split(/; */,$s);
  if ( scalar(@entries) == 1 ) {
    return $entries[0];
  } else {
    return [@entries];
  }
}


sub unparse_multiple {
  my ($x) = @_;
  if (!ref($x)) {
    return $x;
  } else {
    return join(";", @$x);
  }
}

# ------------------------------------------------------------------------

sub parse_indexed {
  my ($s) = @_;
  my @results;
  foreach my $entry ( split(/; */,$s) ) {
    ($entry =~ /^([0-9]+):(.*)/) || die "parse_indexed(... <<$entry>> ...),";
    push @results, [$1,$2];
  }
  return [@results];
}

sub unparse_indexed {
  my ($a) = @_;
  my @in = @$a;
  @in = sort { ($a->[0] <=> $b->[0]) || ($a->[1] cmp $b->[1]) } @in;
  my @out;
  foreach my $entry (@in) {
    my ($pos,$name) = @$entry;
    push @out, $pos.":".$name;
  }
  return join(";",@out);
}

# ------------------------------------------------------------------------

sub parse_named {
  my ($s) = @_;
  my $results = {};
  foreach my $entry ( split(/; */,$s) ) {
    ($entry =~ /^([^:]+):(.*)/) || die "<<$entry>>,";
    $results->{$1} = $2;
  }
  return $results;
}

sub unparse_named {
  my ($h) = @_;
  my @l;
  foreach my $name (sort(keys(%$h))) {
    my $value = $h->{$name};
    push @l, $name.":".$value;
  }
  return join(";",@l);
}

# ------------------------------------------------------------------------
# Translate AA coordinates relative to a specific program to genome
# coordinates relative to the proteins replicon.
# ------------------------------------------------------------------------

sub xlate_aa_coords {
  my ($protein_metadata,$aa_offset,$aa_length) = @_;

  my $protein = $protein_metadata->{protein};
  ( defined($protein_metadata) ) || die "Not found in DB: $protein,";

  my $N_terminal_MAP_processed = $protein_metadata->{N_terminal_MAP_processed};

  my $orf_start = $protein_metadata->{start};
  my $orf_end = $protein_metadata->{end};
  my $orf_strand = $protein_metadata->{strand};

  my $bp_offset = 3 * $aa_offset;
  my $bp_length = 3 * $aa_length;

  my $n_term_fix = ( $N_terminal_MAP_processed ) ? 3 : 0;

  my ($start2,$end2,$strand2);

  if ( $orf_strand eq "+") {
    my $base = $orf_start + $bp_offset + $n_term_fix;
    ($start2,$end2,$strand2) =
      ( $base, $base+$bp_length-1, "+");
  } elsif ( $orf_strand eq "-") {
    my $base = $orf_end - $bp_offset - $n_term_fix;
    ($start2,$end2,$strand2) =
      ( $base-$bp_length+1, $base, "-");
  } else {
    die "orf_strand=<<$orf_strand>>,";
  }

  return ($start2,$end2,$strand2);
}

# ------------------------------------------------------------------------
# Sequence buckets
# ------------------------------------------------------------------------

my $min_subseq_len;

my %sequence_of;
my %buckets;

sub init_buckets {
  my ($n) = @_;
  $min_subseq_len = $n;
}

sub add_sequence {
  my ($id,$seq) = @_;
  $sequence_of{$id} = $seq;
  for (my $i = 0; $i+$min_subseq_len <= length($seq); $i++) {
    my $s = substr($seq,$i,$min_subseq_len);
    my $a = $buckets{$s};
    if (!defined($a)) {
      $buckets{$s} = $a = [];
    }
    push @{$a}, [$id,$i];
  }
}

sub search_sequences {
  my ($subseq) = @_;
  my @results;
  (length($subseq) >= $min_subseq_len) || die "hash size too big,";
  my $s = substr($subseq,0,$min_subseq_len);
  my $a = $buckets{$s};
  if (!defined($a)) {
    return ();
  }
  foreach my $e ( @$a) {
    my ($id,$offset) = @{$e};
    my $sequence = $sequence_of{$id};
    if ($subseq eq substr($sequence,$offset,length($subseq))) {
      push @results, [$id, $offset]
    }
  }
  return @results;
}

use Statistics::Descriptive;

sub dump_bucket_stats {
  print STDERR "# Computing bucket stats\n";

  my @sizes;
  foreach my $a (values(%buckets)) {
    push @sizes, scalar(@$a);
  }

  my $stat = Statistics::Descriptive::Full->new();
  $stat->add_data(@sizes);
  print STDERR "# - num buckets: ", $stat->count(), "\n";
  print STDERR "# - mean bucket size: ", $stat->mean(), "\n";
  print STDERR "# - stdev: ", $stat->standard_deviation(), "\n";
  print STDERR "# - min: ", $stat->quantile(0), "\n";
  print STDERR "# - 25 \%tile: ", $stat->quantile(1), "\n";
  print STDERR "# - 50 \%tile: ", $stat->quantile(2), "\n";
  print STDERR "# - 75 \%tile: ", $stat->quantile(3), "\n";
  print STDERR "# - max: ", $stat->quantile(4), "\n";
}


# ------------------------------------------------------------------------
# comparing psm's
# ------------------------------------------------------------------------

use Data::Dumper;

sub cmp_psms {
  my ($x,$y) = @_;
  # fixme: handle case when ->{spectrum} are array refs
  return ($x->{spectrum} cmp $y->{spectrum})
    || ($x->{peptide} cmp $y->{peptide});
}


# ------------------------------------------------------------------------
# Done.
# ------------------------------------------------------------------------

1;
