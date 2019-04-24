BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

open IN, $ARGV[0];
while (my $l = <IN>) {
  next if ($l =~ /^\!/);

  chomp $l;

  my @a = split /\t/, $l, -1;
  
  my $id = $a[10];

  my @b  = split /\|/,$id; 

  next if ($b[0] eq "");

  print "$b[0]\t$a[4]\n";
  
}
