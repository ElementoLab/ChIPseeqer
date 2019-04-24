#
# run FIRE with 4/5 of the data
# build classifier based on motifs
# estimate performance based on 1/5
#
 

use lib "$ENV{FIREDIR}/SCRIPTS";
use Getopt::Long;
use Sets ;
use Data::Dumper ;

use Table ;
use strict;
use Subs ;
use Bayes ;

my $expfile   = undef ;
my $profiles  = undef ;
my $fastafile = undef;

GetOptions ('expfile=s'     => \$expfile,
	   'fastafile=s'    => \$fastafile);

my $outdir = "$expfile\_CROSS";
mkdir $outdir if (! -e $outdir);


use Table;

my $ta = Table->new;
$ta->loadFile($expfile);
my $a_ref = $ta->getArray();
my $h     = shift @$a_ref;

my $h_ref_p = subdivide_expfile($a_ref, 5);

foreach my $k (keys(%$h_ref_p)) {
  #print "$k\t$h_ref_p->{$k}\n";
}

for (my $i=0; $i<5; $i++) {

  # bulild a new expfile, omiting genes in bin $i  
  open OUT1, ">$outdir/$i.txt";
  print OUT1 join("\t", @$h) . "\n";

  open OUT2, ">$outdir/$i.t.txt";
  print OUT2 join("\t", @$h) . "\n";

  foreach my $r (@$a_ref) {
    #print "$r->[0]\t$h_ref_p->{$r->[0]}\n";
    if ($h_ref_p->{$r->[0]} != $i) {
      print OUT1 "$r->[0]\t$r->[1]\n";
    } else {
      print OUT2 "$r->[0]\t$r->[1]\n";
    }
  }

  close OUT1;
  close OUT2;

  # run FIRE
  my $todo = "fire.pl --expfile=$outdir/$i.txt --fastafile_dna=$fastafile --nodups=1 --minr=2.0";
  system($todo);

  # train Naive Bayes
  $todo = "perl $ENV{FIREDIR}/MORESCRIPTS/naiveBayes.pl --expfile=$outdir/$i.txt --testfile=$outdir/$i.t.txt | columns.pl 1 2 | sort | uniq -c";
  system($todo);

}

#
# 
#
sub subdivide_expfile {

  my ($a_ref, $p) = @_;

  # 1. how many classes ?
  #my %H1 = ();
  my %H2 = ();
  foreach my $r (@$a_ref) {
    #$H1{$r->[1]} ++;
    push @{ $H2{$r->[1]} }, $r->[0];
  }
  my $ne = scalar( keys( %H2 ) );
  my $ng = @$a_ref;

  

  my %OUT = ();

  foreach my $k (keys(%H2)) {
   

 
    my $n = scalar( @{ $H2{$k} } );

    #print "$n GENES with class $k\n";

    my $c = Sets::shuffle_array( $H2{$k} );
    #print scalar(@$c) . " GENES with class $k\n";

    my $b = int($n / $p);

    #my $cnt
    for (my $i=0; $i<$p; $i++) {
      for (my $j=$i*$b; $j<$i*$b+$b; $j++) {
	#print "$c->[$j]\t$k\t$i\n";
	$OUT{ $c->[$j] } = $i;
      }
      if ($i == $p-1) {
	for (my $j=$i*$b+$b; $j<$n; $j++) {
	  
	  #print "$c->[$j]\t$k\t$i\n";
	  $OUT{ $c->[$j] } = $i;
	}
      }
    }

  }
  
  return \%OUT;

}
