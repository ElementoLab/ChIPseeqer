#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Getopt::Long;

if (@ARGV == 0) {
  die "Args: --scanacefile=s --expfile=s --bin=s\n";
}

GetOptions("scanacefile=s"    =>  \$scanacefile,
	   "expfile=s"        =>  \$expfile,
	   "bin=s"            =>  \$bin);


my $ta = Table->new;

# load expfle
$ta->loadFile($expfile);
my $h_ref = $ta->getIndexKV(0,1);

# motif matches
$ta->loadFile($scanacefile);
my $a_ref = $ta->getArray();
foreach my $r (@$a_ref) {
  if (defined($h_ref->{ $r->[0] }) && ($h_ref->{$r->[0]} == $bin)) {
    print join("\t", @$r) . "\n";
  }
}

