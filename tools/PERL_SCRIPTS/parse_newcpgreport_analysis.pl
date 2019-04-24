#!/usr/bin/perl

open IN, $ARGV[0];
while (my $l = <IN>) {
  chomp $l;
  if ($l =~ /^ID\ {3}(.+?)\ /) {
    $id = $1;
  } 

  if ($l =~ /^FT   no islands detected/) {
    #print "$id\t0\n";
  }

  if ($l =~ /FT   CpG island       (\d+)\.\.(\d+)/) {
    
    #print "$id\t$1\t$2\n";
    
    my $tt = "$1-$2";
    push @T, $tt;
    
  }

  if ($l =~ /^\/\//) {
    if (scalar(@T) > 0) {
      print "$id\t1\t"; print join("\t", @T); print "\n";
    } else {
      print "$id\t0\n";
    }
    @T = ();
  }
  

}
close IN;
