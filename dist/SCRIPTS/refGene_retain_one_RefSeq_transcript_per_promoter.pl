#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";

use Sets;

my $refseq = $ARGV[0];

my %ORFPROM = ();
my %ORFNM   = ();

open IN, $refseq;
my %NM = ();
while (my $l = <IN>) {  
  chomp $l;
  my @a = split /\t/, $l, -1;

  next if ($a[2] =~ /\_/);

  my $orf = $a[12];
  my $nm  = $a[1];

  if (defined($ORFNM{$orf}{$nm})) { # we have seen this NM already
    next;
  }
  $ORFNM{$orf}{$nm} = 1;

  my $tss = undef;
  if ($a[3] eq "+") {
    $tss = $a[4];
  } else {
    $tss = $a[5];
  }

  

  push @{ $ORFPROM{$orf}{$tss} }, $nm;
  $NM{$nm} = \@a;

}
close IN;



foreach my $g (keys(%ORFPROM)) {
  
  foreach my $tss (keys(%{$ORFPROM{$g}})) {
    my $small = Sets::get_smallest_NM($ORFPROM{$g}{$tss});
    #print "$g\t$tss\t$small\t" . join("\t", @{$ORFPROM{$g}{$tss}}) . "\n";
    if ($ARGV[1] eq "") {
      print join("\t", @{$NM{$small}}) . "\n";
    } else {
      my @a =  @{$NM{$small}};
      my $st = undef;
      if ($a[3] eq '+') {
	$st = 1;
      } else {
	$st = -1;
      }
      print "$a[1]\t$a[2]\t$a[4]\t$a[5]\t$st\n";
    }
  }

}
