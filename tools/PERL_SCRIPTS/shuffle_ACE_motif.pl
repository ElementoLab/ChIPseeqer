BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;
use strict;




my $perm = undef;

open IN, $ARGV[0];
my $r = <IN>;
print $r;

my $lcnt = 0;
while (my $r = <IN>) {
  chomp $r;

  if ($lcnt == 0) {
    my $m = length($r);    
    $m --;
    $perm  = Sets::shuffle_array([0..$m]);
  }

  my @a = split //,$r,-1; 
  
  for (my $i=0; $i<@$perm; $i++) {
    print $a[$perm->[$i]];
  }
  print "\n";
  $lcnt ++;
} 

close IN;
