#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use strict;

use Getopt::Long;

my $col  = 0;
my $head = 1;



my $ta = Table->new;
$ta->loadFile("$ENV{HOME}/PROGRAMS/ChIPseeqer/DATA/hg18/refLink.txt.07Jun2010.colreag");
my $h_ref = undef;
if ($ARGV[1] eq "") {
  $h_ref = $ta->getIndex(0);
} else {
  $h_ref = $ta->getIndex(1);
}

$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray;

if ($head == 1) {
  my $r = shift @$a_ref;
  print join("\t", @$r);
  if ($ARGV[1] eq "") {
    print "\tDescription\n";
  } else {
    print "\tORF\tDescription\n";
  }
}

foreach my $r (@$a_ref) {
  if ($ARGV[1] eq "") {
    print join("\t", @$r) . "\t" . $h_ref->{$r->[$col]}->[2] . "\n"; 
  } else {
    my $cc = $r->[$col];
    $cc =~ s/\-\d+$//;
    print join("\t", @$r) . "\t" . $h_ref->{$cc}->[0] . "\t" .  $h_ref->{$cc}->[2]  . "\n";
  }
}

