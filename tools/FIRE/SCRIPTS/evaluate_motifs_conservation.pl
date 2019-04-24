use lib "$ENV{FIREDIR}/SCRIPTS";

use Table;
use Sets;

use strict;

if (@ARGV == 0) {
  die "Usage: perl evaluate_motifs_conservation.pl -fasta1 FILE -fasta2 FILE -nbgenes INT -summaryfile FILE -kmerfile FILE -outfile FILE\n";
}

my $progdir     = "$ENV{FIREDIR}/PROGRAMS";

my $fasta1      = Sets::get_parameter(\@ARGV, "-fasta1");
my $fasta2      = Sets::get_parameter(\@ARGV, "-fasta2");
my $nbgenes     = Sets::get_parameter(\@ARGV, "-nbgenes");
my $kmerfile    = Sets::get_parameter(\@ARGV, "-kmerfile");
my $summaryfile = Sets::get_parameter(\@ARGV, "-summaryfile");
my $outfile     = Sets::get_parameter(\@ARGV, "-outfile");
my $shuffle     = undef;
if (Sets::exist_parameter(\@ARGV, "-shuffle")) {
  $shuffle     = Sets::get_parameter(\@ARGV, "-shuffle");
}
my $rna         = undef;
if (Sets::exist_parameter(\@ARGV, "-rna") == 1) {
  $rna         = Sets::get_parameter(\@ARGV, "-rna");
  print "Setting -rna to $rna.\n";
}
my $rootdir     = ".";

#
#  load summary
#
my $ta = Table->new;
$ta->loadFile($summaryfile);
my $a_ref_mo = $ta->getArray();
my %STAT      = ();
my @MOTIFS    = ();
foreach my $r (@$a_ref_mo) {
  push @MOTIFS, $r->[0];
}


#
#  load conservation file
#
$ta->loadFile($kmerfile);
my $a_ref_km = $ta->getArray();

open OUT, ">$outfile" or die "cannot open $outfile\n";
my @HITS = ();
foreach my $re (@MOTIFS) {

  my $c_index = &eval_index($re, $a_ref_km, $fasta1, $fasta2, $nbgenes, $rna);
  my $txt = "$re\t$c_index";
  
  
  if (defined($shuffle) && ($shuffle == 1)) {
    my @a = ();
    my $cnt = 0;
    for (my $i=0; $i<100; $i++) {
      #print "i=$i\n";
      my $myre      = join("", @{ Sets::shuffle_re($re) });
      my $c_index_r = &eval_index($myre, $a_ref_km, $fasta1, $fasta2, $nbgenes, $rna);
      push @a, $c_index_r;
      if ($c_index > $c_index_r) {
	$cnt ++;
      }
    }
    
    my $avg = Sets::average(\@a);
    my $std = Sets::stddev(\@a);
    
    my $z    = sprintf("%3.2f", ($c_index - $avg) / $std);
    my $cntf = sprintf("%3.2f", $cnt/100);
    
    $txt .= "\t$cntf\t$z";
    
  } 

  $txt .= "\n";
  print OUT $txt;
  print $txt;
}

close OUT;
 

sub eval_index {
  my ($re, $a_ref_km, $fasta1, $fasta2, $nbgenes, $rna) = @_;
  
  my $nk       = scalar(@$a_ref_km);
  my $todo = "$progdir/recompare -re $re -fasta1 $fasta1 -fasta2 $fasta2 -nbgenes $nbgenes ";
  if (defined($rna)) {
    $todo .= " -rna $rna ";
  }
  my $out = `$todo`;
  my $cnt = 0;
  if ($out !~ / inf/) {
    my ($lp) = $out =~ / = (.+)$/;
    while (($a_ref_km->[$cnt]->[4] > $lp) || ($a_ref_km->[$cnt]->[4] =~ /inf/)) {
      $cnt ++;
    }
  }

  return sprintf("%3.2f", 100 * (($nk - $cnt) / $nk));
}
