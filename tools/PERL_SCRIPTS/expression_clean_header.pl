#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

open IN, $ARGV[0];

my $l = <IN>;
chomp $l;
my @a = split /\t/, $l, -1;

foreach my $r (@a) {
  $r =~ s/\ .+$//;
}

print join("\t", @a) . "\n";

while (my $l = <IN>) {
  print "$l";
}

close IN;



