#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;
$ta->loadFile("/Users/olivier/PROGRAMS/PLASMODIUM/GENOMES/annotations/pf_annotation_6.0_genedesc.txt");
my $h_ref = undef;
if ($ARGV[1] eq "") {
  $h_ref = $ta->getIndex(0);
} else {
  $h_ref = $ta->getIndex(1);
}

$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray;


foreach my $r (@$a_ref) {
  my $ann  = $h_ref->{$r->[0]}->[2]; 

  if ($ann =~ s/unknown function//) {
    $ann = "UF";
  }
  if ($ARGV[1] eq "") {
    print join("\t", @$r) . "\t" . $ann . "\n"; 
  } else {
    print join("\t", @$r) . "\t" . $h_ref->{$r->[0]}->[0] . "\t" .  $ann  . "\n";
  }
}

