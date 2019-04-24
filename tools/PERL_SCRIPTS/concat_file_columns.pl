#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use FileHandle;

my @FH = ();
foreach my $f (@ARGV) {
  my $fh = IO::File->new($f);
  push @FH, $fh;
}


while (1) {
  my @a = ();
  my $o = 0;
  foreach my $fh (@FH) {
    my $l = <$fh>;
    if (!defined($l)) {
      $o = 1;
      last;
    }
    chomp $l;
    push @a, $l;
  }
  if ($o == 1) {
    last;
  }
  print join("\t", @a); print "\n";
}
  



