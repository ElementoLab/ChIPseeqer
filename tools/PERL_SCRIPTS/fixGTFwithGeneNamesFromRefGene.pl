#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;

my %NM2ORF = ();
open IN, $ARGV[1] or die "Cannot open $ARGV[1]\n";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;

  $NM2ORF{$a[1]} = $a[12];

}
close IN;


# read GTF
open IN, $ARGV[0] or die "Cannot open $ARGV[0]\n";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  
  my $desc = $a[8];

  my ($nm) = $desc =~ /gene_id \"(.+?)\"/;
 
  if (defined($NM2ORF{$nm})) {
    $desc =~ s/gene_id \"(.+?)\"/gene_id \"$NM2ORF{$nm}\"/;
    $desc = "$desc" . "g_id \"$NM2ORF{$nm}\";";
  } else {
    $desc = "$a[8]" . "g_id \"$nm\";";
  }

  $a[8] = $desc;

  print join("\t", @a) . "\n";

}
close IN;

