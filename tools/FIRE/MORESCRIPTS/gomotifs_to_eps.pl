use Table;
use Sets;
use Fire;

use lib qw(PostScript-Simple-0.07/lib);
use PostScript::Simple;

use strict;


if (@ARGV == 0) {
  die "Usage: perl gomotifs_to_eps.pl gofile summary clusters\n";

}
#
# load GO annot
#
my $h_ref_go  = Fire::loadFireGOMotifs($ARGV[0]);

#
# load summary
#
my $h_ref_sum = Fire::loadFireMotifSummary($ARGV[1]);

#
# load clusters to get order
#
my $ta = Table->new;
$ta->loadFile($ARGV[2]);
my $a_ref_clu = $ta->getArray();

#
#  START DRAWING
#
#

my $xbase = 50;
my $ybase = 50;
my $h     = 30;
my $xsize = 595;

my $nm    = @$a_ref_clu;
my $ysize = $nm * $h + $ybase + 10;

my $p = new PostScript::Simple(xsize     => $xsize,
			       ysize     => $ysize,
			       colour    => 1,
			       eps       => 1,
			       units     => "pt");

$p->setlinewidth(0.5);
$p->setcolour("black");
$p->setfont("Times", 12);

$p->text({align => "center"},  $xbase  + 50/2,  $ysize - ( $ybase - 15), "Motif" );

$p->text({align => "center"},  $xbase  + 50 + 35,  $ysize - ( $ybase - 15 ), "Location" );

$p->text({align => "left"},  $xbase + 50 + 70,  $ysize - ( $ybase - 15 ), "Most significant functional enrichment" );

$p->line( $xbase, $ysize - ( $ybase - 10 ), $xbase + 400, $ysize - ( $ybase - 10 ));

my $i = 0;
foreach my $r (@$a_ref_clu) {

  
  my $go  = $h_ref_go->{$r->[0]};

  next if ($go eq "");

  
  my $e1 = &get_eps_logo_object($r->[0], $h_ref_sum->{$r->[0]}->{"RNA"}, 0);
  my $eh = $e1->height;
  my $ew = $e1->width;
  my $ew_new = int(0.5 + $ew * $h / $eh); print "e=$ew_new\n";
  $e1->scale($h / $eh);
  
  $p->_add_eps($e1, $xbase,  $ysize - ( $ybase + $i*$h+$h ) ); 

  my $loc = ($h_ref_sum->{$r->[0]}->{"RNA"}==0?"5'":"3'UTR");

  $p->text({align => "center"},  $xbase  + $ew_new + 35,  $ysize - ( $ybase + $i*$h+$h/2 ), $loc );

  $p->text({align => "left"},  $xbase + $ew_new + 70,  $ysize - ( $ybase + $i*$h+$h/2 ), $go );

  $i ++;
}


$p->output("$ARGV[0].eps");




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
  system("./weblogo/seqlogo -f $d/$cnt.txt -F EPS  -a -c -M -n -Y -w 5 -h 3 > $d/$cnt.eps");

  #
  # START: integrate Motif LOGO
  #  
  my $e  = new PostScript::Simple::EPS(file => "$d/$cnt.eps");
 
  return $e;  
}

