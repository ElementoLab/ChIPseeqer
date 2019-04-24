#!/usr/bin/perl

#
# simple script that summarizes a set of PAGE analyses in a single matrix
#

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;
use strict;

use Getopt::Long;

if (@ARGV == 0) {
  die "Args: --expfiles=PATTERN --draw=INT --tsignf=FLT\n";
}

my $expfiles = undef;
my $draw     = 1;
my $tsignif  = 0.01;
my $outfile  = "combined_matrix";
my $suffix   = undef;

GetOptions("expfiles=s" => \$expfiles,
           "draw=s"     => \$draw,
	   "tsignif=s"  => \$tsignif,
	   "outfile=s"  => \$outfile,
	   "suffix=s"   => \$suffix);
	   

if (defined($suffix)) {
  $outfile .= $suffix;
}

my $a_ref_expfiles = Sets::getFiles($expfiles);

print STDERR "Found " . scalar(@$a_ref_expfiles) . " expfiles\n";

my $ta = Table->new;

my %PV = ();
foreach my $f (@$a_ref_expfiles) {

  my $ff = "$f\_PAGE/pvmatrix.txt";

  $ta->loadFile($ff);
  my $a_ref = $ta->getArray();
  shift @$a_ref;
  foreach my $r (@$a_ref) {
    my $g = shift @$r;
    shift @$r;
    my $p = shift @$r;
    my @a = split /\//, $p;
    my $v = $a[0];
    $PV{ $g }{ $f } = $v;
  }
  
}


open OUT, ">$outfile" or die "Cannot open\n";

print OUT "GO\t" . join("\t", @$a_ref_expfiles) . "\n";


foreach my $g (keys(%PV)) {
  
  

  my $txt = "";
  my $sig = 0;
  foreach my $e (@$a_ref_expfiles) {
    my $p = $PV{$g}{$e};
    if (!defined($p)) {
      #$p = "NA"; 
      $p = sprintf("%4.3f", Sets::log10(0.5));
    } else {

      if ($p < Sets::log10($tsignif)) {
	$sig = 1;
      }
      
      
      #$p = Sets::round_pvalue_up(exp($p * log(10)));

    }
    #my $maxp = 1e-50;
    #if ($p < Sets::log10($maxp)) {
    #  $p = Sets::log10($maxp);
    #}
    $p = -1 * $p;
    $txt .= "\t" . $p;
  }

  if ($sig == 1) {
    print OUT "$g$txt\n";
  }	


}

close OUT;
print "Created $outfile !\n";

if ($draw == 1) {
  my $todo = "draw_overrep_heatmap.pl --coldrawmotifs=0 --matrix=$outfile --draw=open --clustcols=1 --clustrows=1 --drawhorizscale=1 --cmap=/Users/olivier/PROGRAMS/MYTREEVIEW/HEATMAPS/rb_cmap.csv --drawonlysignif=1 --ybase=200 --minmax=0,5 --headerangle=45 --headerfontsize=8 --euclidean=1 --tsignif=$tsignif --algoclust=max";
  print "$todo\n";
  system($todo);
}
