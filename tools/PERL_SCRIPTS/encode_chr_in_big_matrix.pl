BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

open IN, $ARGV[0];
my $n = $ARGV[1];
if ($n eq "") {
  die "Provide col #\n";
}
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  
  $a[$n] =~ s/[Cc]hr//;
  if (uc($a[$n]) eq 'X') {
    $a[$n] = 23;
  }
  if (uc($a[$n]) eq 'Y') {
    $a[$n] = 24;
  }
  
  print join("\t", @a) . "\n";
  

}
close IN;

