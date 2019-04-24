#!/usr/bin/perl

use lib qw(/home/elemento/PERL_MODULES);

use Sets;


#my $a_ref = Sets::readSet($ARGV[0]);
my $h_ref = Sets::getIndex($ARGV[0]);


open IN, $ARGV[1];
if (!defined($ARGV[2])) {
my $l = <IN>;
print "$l";
}
while (my $l = <IN>) {
  chomp $l;  
  my @a = split /\t/, $l, -1;  
  if (!defined($h_ref->{ $a[0] })) { 
    print "$l\n";
  }
      
}
close IN;
