BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use ScanACE;
my $sc = ScanACE->new;


my $a_ref = $sc->_readScanaceResults($ARGV[0]);

my %COUNT = ();
foreach my $r (@$a_ref) {
  push @{ $COUNT{ $r->[0] } },  $ARGV[1] - $r->[1]; 
}


foreach my $k (keys(%COUNT)) {  
  @{ $COUNT{ $k } } = sort { $a <=> $b } @{ $COUNT{ $k } };  
  print "$k\t";
  print scalar(@{ $COUNT{ $k } }); print "\t";
  print join("\t", @{ $COUNT{ $k } }); 
  print "\n";
}
