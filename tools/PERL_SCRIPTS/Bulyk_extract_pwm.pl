#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use Sets;
use strict;

my $outfile = undef;

open IN, $ARGV[0];

#my $l = <IN>;
#print "$l";
#my @a = split /\t/, $l;

#my $outfile = "$a[1].txt";
#

my $num = $ARGV[1];

my $in = 0;
my $cnt = 1;
while (my $l = <IN>) {
  chomp $l;

  if ($l =~ /^(\d+)/) {
    my $txt = "";
    my @a = split /\t/, $l;
    my $file    = Sets::filename($ARGV[0]);
    my $outfile = "PBM_$a[1]_$file";

    #my $outfile = "$a[1]\_$ARGV[0]";

    my $thisnum = $cnt;
    
    $txt .= "$l\n";

    while (($l = <IN>) && ($l !~ /^Probab/)) {
    }
    $l = <IN>;

    for (my $i=0; $i<4; $i++) {
      $l = <IN>;
      $txt .= "$l";
    }
    #last;

    if ($num == $thisnum) {
      open OUT, ">$outfile";
      print OUT $txt;
      print "Created $outfile with motif $num\n";
      close OUT;
      exit;
    }
    

    $cnt++;  # increase motif  count
  }


}
close IN;

