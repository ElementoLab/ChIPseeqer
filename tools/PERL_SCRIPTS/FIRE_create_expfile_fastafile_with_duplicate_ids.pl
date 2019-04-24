#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use strict;
use Getopt::Long;
use Fasta;
use Table;


my $expfile   = undef;
my $fastafile = undef;



GetOptions("expfile=s"   => \$expfile,
	   "fastafile=s" => \$fastafile);


# load fastafile into a hadh table

my $fa = Fasta->new;
$fa->setFile($fastafile);

my %H = ();
while (my $a_ref = $fa->nextSeq()) {
    my ($n, $s) = @$a_ref;
    $H{$n} = $s;
}



my $ta = Table->new;
$ta->loadFile($expfile);
my $a_ref = $ta->getArray();

my $r = shift @$a_ref;

open OUTS, ">$expfile.dups.fasta";
open OUTE, ">$expfile.dups";

print OUTE join("\t", @$r) . "\n";

my %CNTGENES = ();

foreach my $r (@$a_ref) {
  
  next if (!defined($H{$r->[0]}));

  my $n = undef;
  if (defined($CNTGENES{$r->[0]})) {
    $n = $CNTGENES{$r->[0]};
  } else {
    $n = 0;
  }
  
  print OUTE "$r->[0]-$n\t$r->[1]\n";

  print OUTS ">$r->[0]-$n\n$H{$r->[0]}\n\n";


  $CNTGENES{$r->[0]} ++;

}


close OUTE;
close OUTS;
