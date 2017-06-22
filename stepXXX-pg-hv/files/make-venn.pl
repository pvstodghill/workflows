#! /usr/bin/env perl

use strict;
use warnings FATAL => 'all';
use Carp::Always;

# use FindBin;
# use lib "$FindBin::Bin";
# use Xyzzy;

use constant { TRUE => 1, FALSE => 0 };

# ------------------------------------------------------------------------
# Process the command line
# ------------------------------------------------------------------------

use File::Basename;
use Getopt::Std;

our $opt_o = "venn.png";
our $opt_h;

sub usage {
  my $progname = basename($0);
  print STDERR "Usage: $progname [options] sample1.gff sample2.gff ...\n";
  print STDERR "-h - print help\n";
  print STDERR "-o FILE.xxx - output file (xxx: pdf, png) [venn.png]\n";
  exit(@_);
}

my $stat = getopts('ho:');
if (!$stat) {
  usage(1);
}
if ($opt_h) {
  usage();
}

my $opt_o_fmt;
if ( lc($opt_o) =~ /\.png$/ ) {
  $opt_o_fmt = "png";
} elsif ( lc($opt_o) =~ /\.pdf$/ ) {
  $opt_o_fmt = "pdf";
} else {
  die "Unknown format: <<$opt_o>>,";
}

# ------------------------------------------------------------------------
# Read the GFF files; load into tables.
# ------------------------------------------------------------------------

my %rows;
my %cols;

foreach my $file ( @ARGV ) {
  my $colname = $file;
  $colname =~ s/^.*\///;
  $colname =~ s/\..*$//;
  $cols{$colname} = TRUE;
  open(my $fh, "<", $file) || die "<<$file>>,";
  while (<$fh>) {
    if ( /^#/ ) { next; }
    chomp;
    my ($seqname,$source,$feature,$start,$end,
    $score,$strand,$frame,$attributes) = split(/\t/,$_);
    my $rowname = join(",", $seqname,$start,$end,$strand);
    my $row = $rows{$rowname};
    if (!defined($row)) {
      $row = $rows{$rowname} = {};
    }
    $row->{$colname} = TRUE;
  }
  close $fh;
}

# ------------------------------------------------------------------------
# Write the data files for R.
# ------------------------------------------------------------------------

use File::Temp qw/ :POSIX /;

my $r_input_file = tmpnam();

open(my $fh, ">", $r_input_file) || die;

my @colnames = sort(keys(%cols));
print $fh join("\t",@colnames),"\n";

foreach my $rowname ( sort(keys(%rows)) ) {
  my $row = $rows{$rowname};
  my @results;
  foreach my $colname ( @colnames ) {
    my $x = defined($row->{$colname}) ? "1" : "0";
    push @results, $x;
  }
  print $fh join("\t",@results),"\n";
}

close $fh;


# ------------------------------------------------------------------------
# Run R
# ------------------------------------------------------------------------

use Statistics::R;
my $R = Statistics::R->new();
$R || die;

# v---- this is bohhhhhhhhh-gus!
my $extra_lib_dir = $ENV{HOME}."/opt/R";
if ( -e $extra_lib_dir ) {
  print $R->run(qq{.libPaths(new=c("$extra_lib_dir"))}),"\n";
}
# ^---- this is bohhhhhhhhh-gus!

# setup-bioc limma
$R->run(qq{suppressPackageStartupMessages(library(limma))});
$R->run(qq{data <- read.delim("$r_input_file")});
$R->run(qq{a <- vennCounts(data)});
$R->run(qq{$opt_o_fmt("$opt_o")});
$R->run(qq{vennDiagram(a)});
$R->run(qq{dev.off()});

$R->stop();

unlink($r_input_file);

# ------------------------------------------------------------------------
# Done.
# ------------------------------------------------------------------------
