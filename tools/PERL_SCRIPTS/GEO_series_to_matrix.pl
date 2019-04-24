#!/usr/bin/perl

open IN, $ARGV[0];

my $cnt = 0;
while (my $l = <IN>) {
  chomp $l;
  
  $l =~ s/\n//g;
  $l =~ s/\r//g;

  my @a = split /\t/, $l, -1;
  
 if ($l =~ /^\!/) {
   if (($l =~ /^\!Sample_title/) && ($cnt == 0)) {
     $a[0] = "ID_REF";
     foreach my $r (@a) {
       $r =~ s/\"//g;
     }
     print join("\t", @a) . "\n";
     $cnt = 1;
   } 
   
 } elsif (($l !~ /ID_REF/) && ($l ne "")) {
   
   $a[0] =~ s/\"//g;
   my $n = shift @a;
   print "$n";
   foreach my $r (@a) {
     if ($r ne "") {
       printf("\t%4.3f", $r);
     } else {
       printf("\t");
     }
   }
   print "\n";
 }
}
close IN;
