#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

if (@ARGV == 0) {

   die "Args: all_probe_sets clustering_partition\n";

}


# load set of probe sets

my $ta = Table->new;
$ta->loadFile($ARGV[1]);
my $a_ref = $ta->getArray();


# load existing profile

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $h_ref = $ta->getIndexKV(0,1);

my $r = shift @$a_ref;
print "$r->[0]\tEXP\n";
foreach my $r (@$a_ref) {
  
  if (defined($h_ref->{$r->[0]})) {
    print "$r->[0]\t$h_ref->{$r->[0]}\n";
  } else {
    print "$r->[0]\t0\n";
  }

}



