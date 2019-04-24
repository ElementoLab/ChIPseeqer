BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Fasta;
use Table;
use strict;

my $ta = Table->new;
$ta->loadFile($ARGV[1]);
my $h_ref = $ta->getIndexKV(0,1);


my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

my @ALN = ();
my $seq = undef;
my $cnt = 0;
while (my $a_ref = $fa->nextSeq()) {
  my ($n, $s) = @$a_ref;  
  my @a = split //, $s;
  
  if ($cnt == 0) {
    $seq = $s;
  }

  push @ALN, \@a;
  last;
  $cnt ++;
}

$seq =~ s/\-//g;
my @p_seq = split //, $seq;
my $a_ref_pos = Sets::getREMotifPositions($ARGV[2], $seq);

my @a = @{ shift @ALN };

my $i = 0;
my $j = 0;
my $c = 1;
my $bt = undef;
foreach my $r (@a) {
  
  if ($r ne '-') {
    
    if (Sets::in_array($i, @$a_ref_pos)) {
      
      my $co = 0;
      my $st = "";
      
      for (my $k=$i; $k<$i+$ARGV[3]; $k++) {
	$co += $h_ref->{$k+1};
	$st .= $p_seq[$k];
      }
      $co /= $ARGV[3];

      $bt = substr($seq, $i-125, 250+$ARGV[3]);

      #print "RAP1 ($st, $co)";
      #print " $bt ";
      #print "\n";

      print "$i\tM$c\t" . sprintf("%3.2f", $co) . "\n";
      #print ">M$c\n$bt\n\n";
      $c ++;
    }
    
    

    $i++;
  }
  
  $j++;
}    
