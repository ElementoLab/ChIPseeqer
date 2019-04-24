use lib "$ENV{FIREDIR}/SCRIPTS";
use lib "$ENV{FIREDIR}/SCRIPTS/PostScript-Simple-0.07/lib";

use Table;
use Sets;
use Fire;
use PostScript::Simple;
use Getopt::Long;

use strict;

my $summaryfile = undef;

GetOptions ('summaryfile=s'  => \$summaryfile);


# load summary file
my $a_ref_motifs = Fire::loadFireMotifSummaryArray($summaryfile);

print "Loaded " . scalar(@$a_ref_motifs) . " motifs.\n";

# load optim file
my $a_ref_wm     = []; 

foreach my $m (@$a_ref_motifs) {
  
  my $f = "$summaryfile\.optimization_OUT/$m->{SEED}";

  print "Loadinf $f ... ";
  my $a = Fire::loadFireOptimizedWeightMatrices($f);
  print "Done ($a->[0]->{MI}).\n";

  push @$a_ref_wm, $a->[0];

}


my $n1 = @$a_ref_motifs;
my $n2 = @$a_ref_wm;
if ($n1 != $n2) {
  die "Unequal motif numbers ($n1, $n2).\n";
}


#
#  START DRAWING
#
#

my $xbase =  0;
my $ybase = 50;
my $h     = 30;
my $xsize = 400;

my $nm = @$a_ref_motifs;

my $ysize = $nm * $h + $ybase + 10;

my $p = new PostScript::Simple(xsize     => $xsize,
			       ysize     => $ysize,
			       colour    => 1,
			       eps       => 1,
			       units     => "pt");

$p->setlinewidth(0.5);
$p->setcolour("black");
$p->setfont("Times", 8);


for (my $i=0; $i<$nm; $i++) {
  
  my $r1 = $a_ref_motifs->[$i];
  my $r2 = $a_ref_wm->[$i];

  my $re1 = $r1->{MOTIF};
  my $mi1 = $r1->{MI};
  
  my $rna = $r1->{RNA};

  my $e1  = &get_eps_logo_object($re1, $rna, $i . "1");
  
  my $mi2 = $r2->{MI};
  my $e2  = &get_eps_logo_object_WM($r2->{SITES}, $rna, $i . "2");
  
  my $eh1 = $e1->height;
  my $ew1 = $e1->width;

  my $eh2 = $e2->height;
  my $ew2 = $e2->width;

  my $ew_new1 = int(0.5 + $ew1 * $h / $eh1);
  my $ew_new2 = int(0.5 + $ew2 * $h / $eh2);

  # height must be $h, so scale down to $h = k * LO
  $e1->scale($h / $eh1);
  $e2->scale($h / $eh2);
  
  $p->_add_eps($e1, $xsize/2 - $ew_new1 - 0,  $ysize - ( $ybase + $i*$h+$h ) ); 
  $p->_add_eps($e2, $xsize/2           + 0,  $ysize - ( $ybase + $i*$h+$h ) ); 
 
  $p->text({align => "right"}, $xsize/2 - $ew_new1 - 35,  $ysize - ( $ybase + $i*$h+$h/2 ), "MI=$mi1" );
  $p->text({align => "left"}, $xsize/2 + $ew_new2 + 35,  $ysize - ( $ybase + $i*$h+$h/2 ), sprintf("MI=%5.4f", $mi2) );

  #$p->text({align => "right"}, $xsize/2 - $ew_new - 70,  $ysize - ( $ybase + $i*$h+$h/2 ), $go1 );
  #$p->text({align => "left"},  $xsize/2 + $ew_new + 70,  $ysize - ( $ybase + $i*$h+$h/2 ), $go2 );

  #$i ++;
}


$p->output("out.eps");

system("ps2pdf -dEPSCrop -dAutoRotatePages=/None out.eps out.pdf");

sub get_eps_logo_object_WM {
    
  my ($a_ref_sites, $rna, $cnt) = @_;
  
  #
  # make motif logo
  #
  my $mo = join("\n", @$a_ref_sites) . "\n";

  if ($rna == 1) {
    $mo =~ s/T/U/g;
  }
  
  my $d           = "TMP";
  
  if (!defined($cnt)) {
    $cnt         = 0;
  }

  mkdir $d if (! -e $d);   # to be cleaned

  open OUT, ">$d/$cnt.txt" or die "cannot open $d/$cnt.txt\n";
  print OUT $mo;
  close OUT;
  system("$ENV{FIREDIR}/SCRIPTS/weblogo/seqlogo -f $d/$cnt.txt -F EPS  -a -c -M -n -Y -w 7 -h 3 > $d/$cnt.eps");

  #
  # START: integrate Motif LOGO
  #  
  my $e  = new PostScript::Simple::EPS(file => "$d/$cnt.eps");
 
  return $e;  


}

sub get_eps_logo_object {
  
  my ($re, $rna, $cnt) = @_;
  
  #
  # make motif logo
  #
  my $mo = Sets::myre2wm($re);
  if ($rna == 1) {
    $mo =~ s/T/U/g;
  }
  
  my $d           = "TMP";
  
  if (!defined($cnt)) {
    $cnt         = 0;
  }

  mkdir $d if (! -e $d);   # to be cleaned

  open OUT, ">$d/$cnt.txt" or die "cannot open $d/$cnt.txt\n";
  print OUT $mo;
  close OUT;
  system("$ENV{FIREDIR}/SCRIPTS/weblogo/seqlogo -f $d/$cnt.txt -F EPS  -a -c -M -n -Y -w 5 -h 3 > $d/$cnt.eps");

#9 -> 5
#15 -> ?   ? = 15*5/9 

  #
  # START: integrate Motif LOGO
  #  
  my $e  = new PostScript::Simple::EPS(file => "$d/$cnt.eps");
 
  return $e;  
}


sub get_motif_co_orientation {
  my ($re1, $re2) = @_;

  my $d = "TMP";
  
  my $mo1 = Sets::myre2wm($re1);
  open OUT, ">$d/m1.txt" or die "cannot open $d/m1.txt\n";
  print OUT $mo1;
  close OUT;

  my $mo2 = Sets::myre2wm($re2);
  open OUT, ">$d/m2.txt" or die "cannot open $d/m2.txt\n";
  print OUT $mo2;
  close OUT;

  my $score1 = `./alignace2004/CompareACE $d/m1.txt $d/m2.txt -ss`;
  chomp $score1;

  my $re3 = Sets::getComplement($re2);
  my $mo3 = Sets::myre2wm($re3);
  open OUT, ">$d/m3.txt" or die "cannot open $d/m3.txt\n";
  print OUT $mo3;
  close OUT;

  my $score2 = `./alignace2004/CompareACE $d/m1.txt $d/m3.txt -ss`;
  chomp $score2;
  
  if ($score1 > $score2) {
    return 1;
  } else {
    return -1;
  }	

}

