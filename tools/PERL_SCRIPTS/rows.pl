#!/usr/bin/perl

my %H = ();
foreach my $r (@ARGV) {
  $H{$r} = 1;
}

my $i = 0;
while (my $l = <STDIN>) {

  if (defined($H{$i})) {
    print $l;
  } 

  $i++;
}
