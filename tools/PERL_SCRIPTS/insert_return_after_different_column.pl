#!/usr/bin/perl

my $i    = 0;
my $prev = undef;

while (my $l = <STDIN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;

  next if ($l eq "");

  if (defined($prev) && ($a[0] ne $prev)) {
    print "\n";
  }	

  print "$l\n";

  $prev = $a[0];
  
  $i++;
  
}
