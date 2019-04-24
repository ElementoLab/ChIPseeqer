#!/usr/bin/perl

while (my $l = <STDIN>) {

   chomp $l;
   my @a = split /\t/, $l, -1;
   
   if ($a[$ARGV[0]] >= $ARGV[1] ) {
       print join("\t", @a) . "\n";
   }
}
