#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}

use lib "$home/PERL_MODULES";
use Sets;
#use Table;


my $cnt = 1;
#open IN, $ARGV[0] or die "no such file $ARGV[0]\n";
while (my $l = <STDIN>) {

   chomp $l;
	my @a = split /\ /, $l, -1;
	print "$a[0]\t" if ($ARGV[0]); $cnt++;
    print scalar(@a); print "\n";
}
#close IN;

