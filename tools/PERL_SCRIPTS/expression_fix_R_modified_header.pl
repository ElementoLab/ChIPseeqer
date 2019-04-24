BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

open IN, $ARGV[1];
my $lnew = <IN>;
close IN;

open IN, $ARGV[0];
my $lold = <IN>;
print $lnew;
while (my $l = <IN>) {
  print $l;
}
close IN;


