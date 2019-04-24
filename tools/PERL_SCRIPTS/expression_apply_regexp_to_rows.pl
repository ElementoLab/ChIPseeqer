#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();
shift @ARGV;
my $re = shift @ARGV;

my $i = 0;
foreach my $r (@$a_ref) {

  #if (Sets::in_array($i, @ARGV)) {
    $r->[0] =~ s/$re//g;
  #}
  
  print join("\t", @$r) . "\n";
  
  $i++;
}

