#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;
use strict;


my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref_g = $ta->getArray();
shift @$a_ref_g;

my $h_ref = $ta->getIndexShifted();

my $h_ref_ano = undef;

if (defined($ARGV[2])) {
  $ta->loadFile($ARGV[2]);
  $h_ref_ano = $ta->getIndexKV(2, 1);
}

my @a_p = ();
my $p = $h_ref->{$ARGV[1]};
foreach my $s (@$p) {
  if ($s =~ /\-/) {
    push @a_p, 0;
  } else {
    push @a_p, 1;
  }
}
$p = \@a_p;

my %H = ();
foreach my $r (keys(%$h_ref)) {
  #next if ($r eq $ARGV[1]);  
  my $q = $h_ref->{$r};
  my @a_q = ();
  foreach my $s (@$q) {
    if ($s =~ /\-/) {
      push @a_q, 0;
    } else {
      push @a_q, 1;
    }
    
  }
  $q = \@a_q;
  #print join("\t", @a_q) . "\n";
  $H{ $r } = -Sets::euclidean($p, $q);
}

#foreach my $k (keys(%H)) {
#  print "$k\t$H{$k}\n";
#}

my $o_ref = Sets::hash_order(\%H);
foreach my $k (reverse(@$o_ref)) {
  print "$k\t" . sprintf("%4.3f", $H{$k});
  print "\t" . join("\t", @{$h_ref->{$k}});
  if (defined($ARGV[2])) {
    print "\t" . $h_ref_ano->{$k} . "\n";
  }
  print "\n";
}

