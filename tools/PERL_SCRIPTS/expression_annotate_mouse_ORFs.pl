#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;
$ta->loadFile("$ENV{HOME}/PROGRAMS/SNPseeqer/REFDATA/mm9/refGene.txt.25MAY2010.desc");
my $h_ref = undef;
if ($ARGV[1] ne "") {
  $h_ref = $ta->getIndex(0);
} else {
  $h_ref = $ta->getIndex(1);
}

$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray;
my $r = shift @$a_ref;
print join("\t", @$r);
if ($ARGV[1] eq "") {
  print "\tDescription\n";
} else {
  print "\tORF\tDescription\n";
}


foreach my $r (@$a_ref) {
  if ($ARGV[1] ne "") {
    print join("\t", @$r) .  "\t" . $h_ref->{$r->[0]}->[1] . "\t" . $h_ref->{$r->[0]}->[2] . "\n"; 
    # print join("\t", @$r) .  "\t" . $h_ref->{$r->[0]}->[2] . "\n"; 

  } else {
    print join("\t", @$r) . "\t" . $h_ref->{$r->[0]}->[1] . "\t" .  $h_ref->{$r->[0]}->[2]  . "\n";
  }
}

