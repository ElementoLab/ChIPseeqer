#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


my @a = <STDIN>;
my @b = split /\t/, $a[0], -1;
my @c = ();
foreach my $r (@b) {
  push @c, "I";
}
print join("\t", @c) . "\n";

foreach my $r (@a) {
  print $r;
}



