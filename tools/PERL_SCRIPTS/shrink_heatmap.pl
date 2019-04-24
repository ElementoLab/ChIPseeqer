#!/usr/bin/perl


BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

if (@ARGV == 0) {
  die "Usage: perl shrink_heatmap.pl --matrixfile=FILE --clustfile=FILE\n";
}

use Getopt::Long;
use Table;
use Sets;
use strict;

my $matrixfile = undef;
my $clustfile  = undef;
my $fold       = 5;

GetOptions ('matrixfile=s' => \$matrixfile,
	    'fold=s'       => \$fold,
	    'clustfile=s'  => \$clustfile);
	    

my $ta = Table->new;
$ta->loadFile($clustfile);
my $h_ref_c = $ta->getIndexKV(0,1);


my %ROWS = ();

open IN, $matrixfile;
my $l1 = <IN>;

my $cnt = 0;
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  my $n = $a[0];
  my $c = $h_ref_c->{$n};
  push @{ $ROWS{ $c } }, \@a;
  $cnt ++;
}

close IN;


open OUT, ">$matrixfile.fold$fold" or die "Cannot open average matrix file.\n";
open TUO,  ">$clustfile.fold$fold" or die "Cannot open average cluster file (.kgg).\n";
print OUT $l1;
print TUO "CLUSTER\tIDX\n";

my @LL = sort { $a <=> $b } (keys(%ROWS));

#
# traverse clusters, shuffle, print n/fold
#
foreach my $k (@LL) {

  #
  # elements in cluster
  #
  my $n         = scalar(@{ $ROWS{ $k } });
  my $a_ref_shu = Sets::shuffle_array( $ROWS{$k} );

  my $nf        = int( 0.5 + $n / $fold );

  my $a_ref_sma = [];
  
  for (my $i=0; $i<$nf; $i++) {
    my $m = $a_ref_shu->[$i]->[0];
    print TUO "$m\t$k\n";
    print OUT join("\t", @{ $a_ref_shu->[$i] }) . "\n";
  }

}

close OUT;
close TUO;

