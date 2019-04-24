#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;
use strict;


#
# expression
#
my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref_g = $ta->getArray();

# header
my $a_ref_h = shift @$a_ref_g;
shift @$a_ref_h;
my $h_ref_h = {};
for (my $i=0; $i<@$a_ref_h; $i++) {
  $h_ref_h->{$a_ref_h->[$i]} = $i;
}

#
# read phenotype
#
$ta->loadFile($ARGV[1]);
my $a_ref_p = $ta->getArray();
my @a_p     = ();
foreach my $r (@$a_ref_p) {
  $a_p[ $h_ref_h->{$r->[0]} ] = ($r->[1] ne "NA"?$r->[1]:"");
}

# 

print "PHENO\t" . Sets::jointab(\@a_p);

foreach my $r (@$a_ref_g) {
  my $g = shift @$r;
  my $p = Sets::pearson($r,\@a_p);
  #print Sets::jointab(\@a_p);
  #print "$g\t" . Sets::jointab($r);
  print "$g\t$p\n";
}


