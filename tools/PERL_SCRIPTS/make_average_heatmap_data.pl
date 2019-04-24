#!/usr/bin/perl


BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

if (@ARGV == 0) {
  die "Usage: perl make_average_heatmap_data.pl --matrixfile=FILE --clustfile=FILE\n";
}

use Getopt::Long;
use Table;
use strict;

my $matrixfile = undef;
my $clustfile  = undef;

GetOptions ('matrixfile=s' => \$matrixfile,
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
    
    my $n = shift @a;
    my $c = $h_ref_c->{$n};
    
    push @{ $ROWS{ $c } }, \@a;

    $cnt ++;
}

close IN;


# average MATRIX file
open OUT, ">$clustfile.kggavg_data" or die "Cannot open average matrix file.\n";

# reduce CLUSTERING file
open TUO,  ">$clustfile.kggavg_clusters" or die "Cannot open average cluster file (.kgg).\n";
print OUT $l1;
print TUO "CLUSTER\tIDX\n";

my @LL = sort { $a <=> $b } (keys(%ROWS));

foreach my $k (@LL) {
    
    my $n = scalar(@{ $ROWS{ $k } });
    my $l = scalar(@{ $ROWS{ $k }->[0] });
    
    #print "n=$n\n";

    my @exp = ();
    for (my $i=0; $i<$l; $i++) {
      my $cnt = 0;
      my $sum = undef;
      for (my $j=0; $j<$n; $j++) {
	if (($ROWS{ $k }->[$j]->[$i] ne "") && ($ROWS{ $k }->[$j]->[$i] ne 'NULL')) {
	  $sum += $ROWS{ $k }->[$j]->[$i];
	  $cnt += 1;
	}
      }
      if ($cnt > 0) {
	$exp[$i] = sprintf("%4.3f", $sum / $cnt);
      } else {
	$exp[$i] = "";
      }
      
      
    }
    
    print TUO "$k\t$k\n";

    print OUT "$k\t"; 
    print OUT join("\t", @exp); 
    print OUT "\n";

}

close OUT;
close TUO;

