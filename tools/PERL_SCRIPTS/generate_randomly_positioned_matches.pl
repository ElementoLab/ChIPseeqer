#
#  input : the chr lengths
#          the matches
#

use lib qw(/home/olly/PERL_MODULES);
use Sets;
srand;
my $h_ref_chrs = Sets::getIndexKV($ARGV[0], 0,1);

my %H = ();

$p = 0;
# read matches
while (my $l = <STDIN>) {
  chomp $l;
  
  my @a = split /\t/, $l, -1;
  
  my $d = abs($a[1] - $hisp);
  
  if ($d <= 2) {
    $pos = $myp + $d;
    # get a non random position
  } else {
    # get a random one
    $pos = int(rand($h_ref_chrs->{$a[0]}));
  }
  
  #print "$a[1]\t$a[0]\t$pos\n"; 
  
  
  
  push @{ $H{$a[0]} }, $pos; 
  
  $myp  = $pos;
  $hisp = $a[1];
}



foreach my $k (keys(%H)) {
    my @S = sort { $a <=> $b } @{ $H{$k} };
    foreach my $r (@S) {
	print "$k\t$r\n";
    }
}



