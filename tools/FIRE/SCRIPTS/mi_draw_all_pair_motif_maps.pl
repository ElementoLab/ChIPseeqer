use lib "$ENV{FIREDIR}/SCRIPTS";

#
# go thru the mimatrix output, for each coloc pair, draw a motif map
#


use Sets;
use Getopt::Long;
use Table;

use strict;

my $scriptdir     = "$ENV{FIREDIR}/SCRIPTS";
my $mimatrixfile  = undef;
my $profiles      = undef;
my $expfile       = undef;
my $summaryfile   = undef;
my $seqlen        = undef;
my $rootdir       = undef;
my $rightlabel    = undef;
my $leftlabel     = undef;
my $rna           = undef;
my $fullmatrixfile= undef;
my $fastafile     = undef;
my $h             = 3;

if (@ARGV == 0) {
  die "Usage: perl mi_draw_all_pair_motif_maps.pl --mimatrixfile=FILE --profiles=FILE --expfile=FILE --summaryfile=FILE\n";
}

GetOptions ('mimatrixfile=s'   => \$mimatrixfile,
	    'fullmatrixfile=s' => \$fullmatrixfile,
	    'profiles=s'       => \$profiles,
	    'expfile=s'        => \$expfile,
	    'summaryfile=s'    => \$summaryfile,
	    'seqlen=s'         => \$seqlen,
	    'rootdir=s'        => \$rootdir,
	    'leftlabel=s'      => \$leftlabel,
	    'rightlabel=s'     => \$rightlabel,
	    'fastafile=s'      => \$fastafile,
	    'h=s'              => \$h,
	    'rna=s'            => \$rna
	    );



die "Please define -rna.\n" if (!defined($rna));

my $ta = Table->new;


#
# draw individual motif maps
#

my $d = "$summaryfile\_OUT"; mkdir $d if (! -e $d);

$ta->loadFile($summaryfile);
my $a_ref_su   = $ta->getArray();
my $h_ref_rna  = $ta->getIndexKV(0,1);
my $h_ref_seed = $ta->getIndexKV(0,8);

foreach my $r (@$a_ref_su) {
  
  next if ($r->[1] != $rna);

  my $f = $r->[8];
  if (($r->[8] eq '') || ($r->[8] eq '0')) {
    $f = $r->[0];
  }
  print "Draw motif map for $r->[0] .. ";
  my $cmd = "perl $scriptdir/mi_draw_motif_map.pl --expfile=$expfile --clusters=-1 --profiles=$profiles --rna=$r->[1] --motifs_m=$r->[0] --summaryfile=$summaryfile --overrep=0 --outeps=$d/$f.eps --rightlabel=$rightlabel --leftlabel=$leftlabel --fullmatrixfile=$fullmatrixfile --h=$h";

  if (defined($fastafile)) {
    $cmd .= " --fastafile=$fastafile";
  } else {
    $cmd .= " --seqlen=$seqlen";
  }

  system("$cmd > /dev/null");
  print "Done.\n";
  system ("ps2pdf -dEPSCrop -dAutoRotatePages=/None $d/$r->[8].eps $d/$r->[8].pdf");
}



exit if (!defined($mimatrixfile));

$d = "$mimatrixfile\_OUT"; mkdir $d if (! -e $d);


#
# draw interaction motif maps
#
$ta->loadFile($mimatrixfile);
$ta->sortbycol(4, 0);  # sort by z-score, descending
my $a_ref_mi = $ta->getArray();


my $cnt = 0;

foreach my $p (@$a_ref_mi) {


  next if ($h_ref_rna->{$p->[0]} != $h_ref_rna->{$p->[1]});  # must be of the same type
  next if ($h_ref_rna->{$p->[0]} != $rna); # must respect the type set by $rna
  next if ($p->[7] eq 'nan');
  next if ($p->[6]  > 100);

  print "Draw motif map for $p->[0],$p->[1] .. ";

  my $f = "";
  if (($h_ref_seed->{$p->[0]} ne "") && ($h_ref_seed->{$p->[1]} ne "") && ($h_ref_seed->{$p->[0]} ne "0") && ($h_ref_seed->{$p->[1]} ne "0")) {
    $f = $h_ref_seed->{$p->[0]} . "_" . $h_ref_seed->{$p->[1]};
  } else {
    $f = "$p->[0]_$p->[1]";
  }

  my $cmd = "perl $scriptdir/mi_draw_motif_map.pl --expfile=$expfile --clusters=-1 --profiles=$profiles --seqlen=$seqlen --rna=$rna --motifs_m=$p->[0],$p->[1] --summaryfile=$summaryfile --overrep=0 --intersection=1 --outeps=$d/$f.eps --rightlabel=$rightlabel --leftlabel=$leftlabel --fullmatrixfile=$fullmatrixfile --h=3";
  
  system("$cmd > /dev/null");
 

  #print "ps2pdf -dEPSCrop -dAutoRotatePages=/None $d/$cnt.eps $d/$cnt.pdf\n";
  system ("ps2pdf -dEPSCrop -dAutoRotatePages=/None $d/$f.eps $d/$f.pdf");
  
  print " ($d/$f.eps) Done.\n";
  

  
  $cnt ++;
}





