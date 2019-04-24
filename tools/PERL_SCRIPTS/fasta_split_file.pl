#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Fasta;
use strict;

die "Args: fasta num\n" if (@ARGV == 0);

my $fafa = $ARGV[0];

# count numreads
open IN, $fafa or die "Cannot open $fafa\n";
my $numreads = 0;
while (my $l = <IN>) {
  $numreads++ if ($l =~ /^\>/);
}
close IN;

print STDERR "Found $numreads reads.\n";

die "Please provide numfiles\n" if ($ARGV[1] eq "");

my $numreadsperfile = int(0.5 + $numreads / $ARGV[1] );

my $fa = Fasta->new;
$fa->setFile($fafa);

my $cnt = 0;
my $id  = 0;
while (my $a_ref = $fa->nextSeq()) {
  my ($n, $s) = @$a_ref;
  if ($cnt % $numreadsperfile == 0) {
    close OUT if ($cnt > 0);      
    open OUT, ">$fafa.$id";
    print STDERR "$cnt reads read, creating $fafa.$id\n";
    $id ++;
  }
  print OUT ">$n\n$s\n";
  $cnt ++;
}
close OUT;
