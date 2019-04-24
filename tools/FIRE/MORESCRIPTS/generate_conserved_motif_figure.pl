use Table;
use Sets;
use Fire;

use lib qw(PostScript-Simple-0.07/lib);
use PostScript::Simple;

use strict;

# load list of spe1 / spe2 motifs
my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref_m = $ta->getArray();

# load names 1/2 
my $h_ref_names1 = Fire::loadFireMotifNames($ARGV[1]);
my $h_ref_names2 = Fire::loadFireMotifNames($ARGV[2]);

# load GO 1/2
my $h_ref_go1    = Fire::loadFireGOMotifs  ($ARGV[3]);
my $h_ref_go2    = Fire::loadFireGOMotifs  ($ARGV[4]);

# load restrict

if (defined($ARGV[5])) {
  my $a_ref_rest   = Sets::readSet($ARGV[5]);
  my %H = (); my $i = 0;
  foreach my $r (@$a_ref_m) {
    $H{$r->[0]} = $i++;
  }
  
  my @a = ();
  foreach my $r (@$a_ref_rest) {
    push @a, $a_ref_m->[ $H{$r} ];
  }
  @$a_ref_m = @a;
}


#
#  START DRAWING
#
#

my $xbase =  0;
my $ybase = 50;
my $h     = 30;
my $xsize = 795;

my $nm    = @$a_ref_m;
my $ysize = $nm * $h + $ybase + 10;

my $p = new PostScript::Simple(xsize     => $xsize,
			       ysize     => $ysize,
			       colour    => 1,
			       eps       => 1,
			       units     => "pt");

$p->setlinewidth(0.5);
$p->setcolour("black");
$p->setfont("Times", 12);

$p->text({align => "center"}, $xsize/2 - 50/2,  $ysize - ( $ybase - 20  ), "Human" );
$p->text({align => "center"}, $xsize/2 - 50/2,  $ysize - ( $ybase -  8  ), "motif" );



$p->text({align => "center"}, $xsize/2 + 50/2,  $ysize - ( $ybase - 20  ), "Mouse" );
$p->text({align => "center"}, $xsize/2 + 50/2,  $ysize - ( $ybase -  8  ), "motif" );

$p->text({align => "center"}, $xsize/2 - 50 - 35,  $ysize - ( $ybase - 10  ), "Name" );
$p->text({align => "center"}, $xsize/2 + 50 + 35,  $ysize - ( $ybase - 10  ), "Name" );

$p->text({align => "right"}, $xsize/2 - 50 - 70,  $ysize - ( $ybase -10 ), "Most signicant functional enrichment" );
$p->text({align => "left"},  $xsize/2 + 50 + 70,  $ysize - ( $ybase -10 ), "Most signicant functional enrichment" );

$p->line($xbase + 20, $ysize - ($ybase - 3) , $xsize - 20, $ysize - ($ybase - 3));

my $i = 0;
foreach my $r (@$a_ref_m) {
  
  my $re1 = $r->[0];
  my $re2 = $r->[1];

  my $si = get_motif_co_orientation($re1, $re2);
  if ($si == -1) {
    $re2 = Sets::getComplement($re2);
  }
  

  my $e1 = &get_eps_logo_object($re1, 0, 0);
  my $e2 = &get_eps_logo_object($re2, 0, 1);
  
  my $eh = $e1->height;
  my $ew = $e1->width;
  my $ew_new = int(0.5 + $ew * $h / $eh);

  # height must be $h, so scale down to $h = k * LO
  $e1->scale($h / $eh);
  $e2->scale($h / $eh);
  
  $p->_add_eps($e1, $xsize/2 - $ew_new - 0,  $ysize - ( $ybase + $i*$h+$h ) ); 
  $p->_add_eps($e2, $xsize/2           + 0,  $ysize - ( $ybase + $i*$h+$h ) ); 
 
  my $na1 = $h_ref_names1->{$r->[0]};
  $na1 =~ s/^M\d{5}\_//;
  $na1 =~ s/^J\_M.\d+\_//;
  $na1 =~ s/\.txt$//;

  my $na2 = $h_ref_names2->{$r->[1]};
  $na2 =~ s/^M\d{5}\_//;
  $na2 =~ s/^J\_M.\d+\_//;
  $na2 =~ s/\.txt$//;

  my $go1 = $h_ref_go1->{$r->[0]};
  my $go2 = $h_ref_go2->{$r->[1]};

  $p->text({align => "center"}, $xsize/2 - $ew_new - 35,  $ysize - ( $ybase + $i*$h+$h/2 ), $na1 );
  $p->text({align => "center"}, $xsize/2 + $ew_new + 35,  $ysize - ( $ybase + $i*$h+$h/2 ), $na2 );

  $p->text({align => "right"}, $xsize/2 - $ew_new - 70,  $ysize - ( $ybase + $i*$h+$h/2 ), $go1 );
  $p->text({align => "left"},  $xsize/2 + $ew_new + 70,  $ysize - ( $ybase + $i*$h+$h/2 ), $go2 );

  $i ++;
}


$p->output("out.eps");




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

