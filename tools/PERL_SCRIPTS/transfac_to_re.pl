BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;

open IN, $ARGV[0];

my @W_RAW = ();
my @nt = ('A', 'C', 'G', 'T');

my @s  = ();

while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  
  if ($l =~ /^P0/) {
    
    my $cntp = 0;
    while (my $m = <IN>) {
      chomp $m;
      last if ($m =~ /^XX/);
      
      my @b = split /\ +/, $m, -1;
      shift @b;
      my $u = pop @b;
      push @s, $u;
      # print "$u\n";

    }
  }



}
close IN;


my $f = $ARGV[0]; $f =~ s/\.txt//;
$f =~ s/\_+$//;
print "$f\t";

my $t   = join("", @s);
my $tre = Sets::consensus2re($t);

print "$tre\t$t\n";
  

