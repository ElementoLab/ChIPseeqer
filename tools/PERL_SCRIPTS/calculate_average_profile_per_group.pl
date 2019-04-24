#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
#use lib "$ENV{PERLMODULESDIR}";
use lib "$ENV{HOME}/PERL_MODULES";

use Table;
use strict;
if (@ARGV == 0) {
  die "Args: groupfile profiles\n";
}
my $ta = Table->new;

# LOAD PROFILES
$ta->loadFile($ARGV[1]);
$ta->processHeader();
my $head  = $ta->getHeader();
my $h_ref = $ta->getIndexShifted();



#
# load groups
#
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my @groups = ();
foreach my $r (@$a_ref) {

  if (defined($h_ref->{$r->[0]})) {
    push @{ $groups[ $r->[1] ] }, $h_ref->{$r->[0]};
  }

}


my $m = @{$groups[0]->[0]};
my @col = ();
for (my $j=0; $j<$m; $j++) {
  push @col, $j;
}
print join("\t", @$head) . "\n";


for (my $k=0; $k<@groups; $k++) {
  

  my $g = $groups[$k];

  next if (!defined($g));
  my $n = @$g;

  my $m = @{$g->[0]};
  
  my @sum_g = ();

  for (my $i=0; $i<$n; $i++) {
    for (my $j=0; $j<$m; $j++) {
      $sum_g[$j] += $g->[$i]->[$j];
    }
  }
    
  for (my $j=0; $j<$m; $j++) {
    my $prec = $ENV{PREC};
    if ($prec eq "") {
       $prec=2;
    }
    $sum_g[$j] = sprintf("%3.$prec" . "f", $sum_g[$j]/$n );
  }

  print "$k\t" . join("\t", @sum_g) . "\n";

  
}


