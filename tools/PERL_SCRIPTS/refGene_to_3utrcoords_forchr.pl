BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

foreach my $r (@$a_ref) {
  next if ($r->[1] ne $ARGV[1]);
  my $st = $r->[4];
  my $s  = undef;
  my $e  = undef;
  if ($st == 1) {    
    $s = $r->[3];  # end protein
    $e = $r->[6];  # end transcript
    if ($e - $s > 5000) {
      $e = $s + 5000;
    }
  } else {    
    $s = $r->[5];  # start transcript on chr  
    $e = $r->[2];  # start protein 
    if ($e - $s > 5000) {
      $s = $e - 5000;
      if ($s < 0) {
	$s = 0;
      }	
    }
  }
  print "$r->[0]\t$s\t$e\t$st\n";
}

