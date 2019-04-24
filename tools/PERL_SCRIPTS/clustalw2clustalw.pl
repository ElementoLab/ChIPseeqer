#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

open IN, $ARGV[0];
my $l = <IN>;
print $l;
while (my $l = <IN>) {

  chomp $l;
  if ($l eq "") {
    print "$l\n";
  } else {	
    my ($n, $s) = $l =~ /^(.+?)\ +(.+)$/;
    print $n . (" " x (16-length($n))) . "$s\n";
  }

}
close IN;
