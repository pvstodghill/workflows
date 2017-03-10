#! /usr/bin/env perl

use strict;
use warnings;

use FindBin;


my %accession_of;
my $scaffolds_file = $FindBin::Bin."/scaffolds.txt";
open(my $fh, "<", $scaffolds_file) || die "Cannot open <<$scaffolds_file>>,";
while ( <$fh> ) {
  chomp;
  my ($accession,$scaffoldId,$empty) = split(" ");
  (defined($scaffoldId) && !defined($empty)) ||
    die "Ill-formed scaffold line, <<$_>>,";
  if ( $accession eq "accession" ) { next; }
  $accession_of{$scaffoldId} = $accession;
}
close $fh;

# ------------------------------------------------------------------------


my ($genes_file,$named_file) = @ARGV;
if (!defined($named_file)) {
  print STDERR "Usage $0 foo.genes foo.named > foo.gff\n";
  exit(1);
}


# ------------------------------------------------------------------------

my %loc;

open(F, "<", $genes_file) || die "Cannot open <<$genes_file>>,";
while (<F>) {
  chomp;
  my ($locusId,$accession_Protein,$GI,$scaffoldId,$start,$stop,$strand,
      $sysName,$name,$desc,$COG,$COGFun,$COGDesc,$TIGRFam,
      $TIGRRoles,$GO,$EC,$ECDesc) = split(/\t/);
  if ( $locusId eq "locusId" ) { next; }
  my $accession = $accession_of{$scaffoldId};
  (defined($accession))
    || die "Update <<$scaffolds_file>> for scaffoldId=<<$scaffoldId>>,";
  my ($begin,$end);
  ($begin,$end) = ($start,$stop) if ( $strand eq "+");
  ($begin,$end) = ($stop,$start) if ( $strand eq "-");
  (defined($begin)) || die "strand=<<$strand>>,";
  ($begin <= $end) || die;
  $sysName =~ s/_//g;
  $loc{$sysName} = [$accession,$begin,$end,$strand];
}
close F;

# ------------------------------------------------------------------------

my %upstream; # ($upstream{A} == B) === B upstream of A
my %downstream;

open(F, "<", $named_file) || die "Cannot open <<$named_file>>,";
while (<F>) {
  chomp;
  my ($Gene1,$Gene2,$SysName1,$SysName2,$Name1,$Name2,
      $bOp,$pOp,$Sep,$MOGScore,$GOScore,$COGSim,
      $ExprSim) = split(/\t/);
  if ( $Gene1 eq "Gene1" ) {
    next;
  }
  if ( $bOp ne "TRUE" ) { next; }
  my $sign = compare_position($SysName1,$SysName2);
  if ( !defined($sign) ) {
    print STDERR "# Skipping edge: $SysName1 <-> $SysName2\n";
    next;
  } elsif ( $sign < 0 ) {
    # $SysName1 upstream of $SysName2
    !(defined($upstream{$SysName2})) || die "SysName2 = $SysName2,";
    $upstream{$SysName2} = $SysName1;
    !(defined($downstream{$SysName1})) || die "SysName1 = $SysName1,";
    $downstream{$SysName1} = $SysName2;
  } elsif ( $sign > 0 ) {
    # $SysName2 upstream of $SysName1
    !(defined($upstream{$SysName1})) || die "SysName1 = $SysName1,";
    $upstream{$SysName1} = $SysName2;
    !(defined($downstream{$SysName2})) || die "SysName2 = $SysName2,";
    $downstream{$SysName2} = $SysName1;
  } else {
    die;
  }
}
close F;

sub compare_position {
  my ($a,$b) = @_;
  $a =~ s/_//g;
  $b =~ s/_//g;
  my $aa = $loc{$a};
  my $bb = $loc{$b};
  if (!defined($aa)) { return; }
  if (!defined($bb)) { return; }
  my ($aa_accession,$aa_begin,$aa_end,$aa_strand) = @$aa;
  my ($bb_accession,$bb_begin,$bb_end,$bb_strand) = @$bb;
  ($aa_accession eq $bb_accession) || die;
  if ( $aa_strand eq "+" && $bb_strand eq "+" ) {
    return $aa_begin <=> $bb_begin;
  } elsif ( $aa_strand eq "-" && $bb_strand eq "-" ) {
    return $bb_begin <=> $aa_begin;
  } else {
    die;
  }
}

# ------------------------------------------------------------------------

my %first; # most upstream gene of an operon
my %is_first;

foreach my $SysNameA ( keys %loc ) {
  compute_first($SysNameA);
  (defined($first{$SysNameA})) || die "SysNameA = $SysNameA,";
}

sub compute_first {
  my ($SysNameA) = @_;
  if ( defined($first{$SysNameA}) ) {
    return;
  }
  my $SysNameB = $upstream{$SysNameA};
  if ( defined($SysNameB) ) {
    compute_first($SysNameB);
    $first{$SysNameA} = $first{$SysNameB};
  } else {
    $first{$SysNameA} = $SysNameA;
    $is_first{$SysNameA} = 1;
  }
}

# ------------------------------------------------------------------------

# foreach my $SysName ( keys %is_first ) {
#   my @l;
#   while ( defined($SysName) ) {
#     push @l, $SysName;
#     $SysName = $downstream{$SysName};
#   }
#   print join(" -> ", @l),"\n";
# }

# ------------------------------------------------------------------------

my $source = "VIMSS";
my $feature = "misc_feature";
my $score = ".";
my $frame = ".";

foreach my $SysName ( keys %loc ) {
  my ($accession,$begin,$end,$strand) = @{$loc{$SysName}};
  my $operon = $first{$SysName};
  my $first = ($SysName eq $operon) ? "TRUE" : "FALSE";
  my $name = $SysName;
  $name =~ s/PSPTO/PSPTO_/;
  $name =~ s/__+/_/g;
  my $operon_name = $operon;
  $operon_name =~ s/PSPTO/PSPTO_/;
  $operon_name =~ s/__+/_/g;
  my $attributes = sprintf "name %s;vimss_operon %s;vimss_first %s", $name, $operon_name, $first;
  print join("\t",$accession,$source,$feature,$begin,$end,
    $score,$strand,$frame,$attributes),"\n";
}

