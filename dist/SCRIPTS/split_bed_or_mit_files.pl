#!/usr/bin/perl

use lib "$ENV{CHIPSEEQERDIR}";

if ($ENV{CHIPSEEQERDIR} eq "") {
   die "Please set CHIPSEEQERDIR variable.\n";
}



use Sets;
use FileHandle;
use strict;



if (@ARGV == 0) {
  die "Args: readfiles [ can be *.bed ]\n";
}

my $verbose = 1;

my %FH = ();
my $num = 0;

foreach my $f (@ARGV) {
  
  print STDERR "Opening $f\n";

  my $dirname = Sets::dirname($f);
  if ($dirname eq "") {
    $dirname = ".";
  }
  print STDERR "Current directory = $dirname\n";

  open IN, $f or die "Cannot open $f\n";
  
  while (my $l = <IN>) {
    next if ($l =~ /^track/);
    chomp $l;

    #print "$l\n" if ($verbose == 1);

    my @a = split /[\ \t]/, $l, -1;
    $l    = join("\t", @a);    
    my $chr = $a[0];
    
    my $fh = undef;
    
    if (!defined($FH{$chr})) {
      print STDERR "Creating $dirname/reads.$chr\n";
      $FH{$chr} = new IO::File ">$dirname/reads.$chr";
    }
    
    $fh = $FH{$chr};
    
    print $fh "$l\n"; 
    
    $num++;
  }
  close IN;

}

print "All done. Type 'ls' to view the files.\n";

