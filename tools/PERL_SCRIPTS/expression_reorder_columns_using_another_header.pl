#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;
use Table;
use strict;
if (@ARGV == 0) {
  die "Args: matrix refmatrix\n";
}
my $ta = Table->new;
$ta->loadFile($ARGV[0]) or "Cannot load matrix1 file ... \n";
my $a_ref1 = $ta->getArray();

open IN, $ARGV[1] or die "Cannot find second file ...\n";
my $l = <IN>;
chomp $l;
my @a = split /\t/, $l, -1;
my %H = ();
#shift @a;
for (my $i=1; $i<@a; $i++) {
  my ($n) = $a[$i]; # =~ /^(\d+)\ \-/;
  $H{ $n } = $i;
  #
}
close IN;


# make order

my @a = @{ $a_ref1->[0] };
my @o = ();
push @o, 0;
#my $n = shift @a;
for (my $i=1; $i<@a; $i++) {
  my ($p) = $a[$i];  # =~ /^(\d+)\ \-/;
  die "Cannot find match for $p ... \n" if (!defined($H{$p}));
  #print "Pushing $H{$p}.\n";
  $o[ $i ] = $H{$p};
  print STDERR "$i => $H{$p} ($p)\n";
}

# change order
foreach my $r (@$a_ref1) {
  my @out = ();
  for (my $i=0; $i<@$r; $i++) {
    $out[ $o[$i] ] = $r->[$i];

  }

  print join("\t", @out) . "\n";
  #exit;
}




