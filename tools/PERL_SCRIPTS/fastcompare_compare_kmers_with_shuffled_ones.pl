BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;
use FileHandle;
use strict;

STDOUT->autoflush(1);

my $ta    = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref      = $ta->getArray();
my $h_ref      = $ta->getIndexKV(0,1);
my $cnt_passed = 0;
my $cnt = 1;
foreach my $r (@$a_ref) {
  my $k       = length($r->[0]);

  my $a_ref_p = Sets::get_all_permutations($k);

  my @a = split //, $r->[0];
  
  my %H = ();
  foreach my $s (@$a_ref_p) {
    my $kmer = "";
    foreach my $t (@$s) {
      $kmer .= $a[$t];
    }
    if (defined($h_ref->{$kmer})) {
      $H{ $kmer } = 1;
    } else {
      #if (!defined($ARGV[1])) { 
      #print "Adding revcomp\n";
      $H{ Sets::getComplement($kmer) } = 1;
      #}
    }
  }


  my @vals = ();
  foreach my $km (keys(%H)) {
    next if ($km eq $r->[0]);
    #print " $km " . $h_ref->{$km} . "\n";
    push @vals, $h_ref->{$km};
  }
  
  my $z = undef;
  if (@vals > 0) {
    
    my $avg = Sets::average(\@vals); #sprintf("%4.3f", 
    my $std = Sets::stddev(\@vals); 
    
   
    
    if ($std < 1e-20) {
      $z = 0;
    } else {
      $z   = sprintf("%4.3f", ($r->[1] - $avg)/$std);
    }

  } else {
    $z = 0;
  }

  my $np = @vals;

  print "$cnt\t$cnt_passed\t$r->[0]\t$z\t$np\n";
    
  if ($z < 3.0) {
    $cnt_passed ++;
  } else {
    $cnt_passed = 0;
  }
  $cnt ++;
  #last if ($cnt_passed == 10);
  
}
