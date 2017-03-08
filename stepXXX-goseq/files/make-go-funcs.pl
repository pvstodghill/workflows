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
  print STDERR "Usage: $progname [options] NC_004578.gff /path/to/go\n";
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

my ($gff_file,$go_dir) = @ARGV;

# ------------------------------------------------------------------------
# Read the feature table.
# ------------------------------------------------------------------------

my %locus_of_geneid; # $geneid -> $locus_tag
my %unclassified; # $geneid -> $unclassified_p

my $fh;
open($fh,"<",$gff_file) || die;
while (<$fh>) {
  if ( /^#/) { next; }
  chomp;
  my ($seqname,$source,$feature,$start,$end,
    $score,$strand,$frame,$attributes) = split(/\t/,$_);
  if ($feature eq "gene") {
    ( $attributes =~ /(^|;)ID=([^;]+)/ ) || die;
    my $gene = $2;
    ( $attributes =~ /(^|;)locus_tag=([^;]+)/ ) || die;
    my $locus = $2;
    ( $attributes =~ /(^|;)Dbxref=GeneID:([^;]+)/ ) || die;
    my $geneid = $2;
    $locus_of_geneid{$geneid} = $locus;
    $unclassified{$geneid} = TRUE;
  }
}
close $fh;

print STDERR scalar(keys(%locus_of_geneid))," gene's\n";

# ------------------------------------------------------------------------
# Read the GO terms for each locus
# ------------------------------------------------------------------------

my %go_of_geneid; # $geneid -> @go_terms

if ( -e "$go_dir/partial-idmapping_selected.tab" ) {
  my $file = "$go_dir/partial-idmapping_selected.tab";
  open($fh,"<",$file) || die "Failed to open file <<$file>>,";
} elsif ( -e "$go_dir/idmapping_selected.tab.gz" ) {
  my $file = "$go_dir/idmapping_selected.tab.gz";
  open($fh, "gunzip -c $file |") || die "Failed to open file <<$file>>,";
} else {
  die "Cannot find any variation of idmapping_selected.tab file,";
}

while (<$fh>) {
  chomp;
  my ($UniProtKB_AC, $UniProtKB_ID, $GeneID, $RefSeq, $GI, $PDB, $GO,
      $f8,$f9,$f10,$f11,$f12,$f13,$f14,$f15,$f16,$f17,$Genbank)
    = split(/\t/);
  if ( $GO eq "" ) { next; }
  foreach my $geneid ( split(/; /,$GeneID) ) {
    ( $geneid =~ /^[0-9]+$/ ) || die "geneid=$geneid,<<$_>>,";
    if ( !defined($locus_of_geneid{$geneid}) ) {
      next;
    }
    $go_of_geneid{$geneid} = [split(/; /,$GO)];
    delete $unclassified{$geneid};
  }
}
close $fh;

print STDERR scalar(keys(%unclassified))," unclassified gene's\n";

foreach my $geneid (keys(%unclassified)) {
  $go_of_geneid{$geneid} = [];
}

# ------------------------------------------------------------------------
# Read the GO nodes and edges
# ------------------------------------------------------------------------

my %go_of_id; # db id -> GO:XXX 
my %id_of_go; # GO:XXX -> db id
my %descr_of_go; # GO:XXX -> human readable term name

# http://geneontology.org/page/lead-database-schema#go-graph.table.term
open($fh,"<","$go_dir/term.txt") || die;
while (<$fh>) {
  chomp;
  my ($id,$name,$type_type,$acc,$is_obsolete,$is_root,$is_relation) = split("\t");
  (!defined($go_of_id{$id})) || die;
  $go_of_id{$id} = $acc;
  (!defined($id_of_go{$acc})) || die;
  $id_of_go{$acc} = $id;
  $descr_of_go{$acc} = $name;
}
close $fh;

my %parents_of_go; # GO:XXX -> parents of GO:XXX
my %children_of_go; # GO:XXX -> children of GO:XXX

# http://geneontology.org/page/lead-database-schema#go-graph.table.term2term
open($fh,"<","$go_dir/term2term.txt") || die;
while (<$fh>) {
  chomp;
  my ($id,$relationship_type_id,$term1_id,$term2_id,$complete) = split("\t");
  my $term1_acc = $go_of_id{$term1_id};
  (defined($term1_acc)) || die;
  my $term2_acc = $go_of_id{$term2_id};
  (defined($term2_acc)) || die;
  if (!defined($parents_of_go{$term2_acc})) {
    $parents_of_go{$term2_acc} = [];
  }
  push @{$parents_of_go{$term2_acc}}, $term1_acc;
  if (!defined($children_of_go{$term1_acc})) {
    $children_of_go{$term1_acc} = [];
  }
  push @{$children_of_go{$term1_acc}}, $term2_acc;
}
close $fh;

# ------------------------------------------------------------------------
# The list of accessors of a GO:XXX accession
# ------------------------------------------------------------------------

my %ancestors_of_go;
sub ancestors_of_go {
  my ($acc) = @_;
  my $h = $ancestors_of_go{$acc};
  if (defined($h)) {
    return $h;
  }
  $h = ancestors_of_go_list(@{$parents_of_go{$acc}});
  $h->{$acc} = TRUE;		# a node is its own ancestor
  $ancestors_of_go{$acc} = $h;
  return $h;
}

sub ancestors_of_go_list {
  my @accs = @_;
  my $h = {};
  foreach my $go ( @accs ) {
    my $h2 = ancestors_of_go($go);
    foreach my $go2 ( keys %$h2 ) {
      $h->{$go2} = TRUE;
    }
  }
  return wantarray ? (keys %$h) : $h;
}

# ------------------------------------------------------------------------
#
# ------------------------------------------------------------------------

foreach my $geneid (sort(keys(%go_of_geneid))) {
  my $locus = $locus_of_geneid{$geneid};
  my @immediate_terms = @{$go_of_geneid{$geneid}};
  my @all_terms = ancestors_of_go_list(@immediate_terms);
  foreach my $go_term ( sort(@all_terms) ) {
    print $locus,"\t",$go_term,"\n";
  }
}
