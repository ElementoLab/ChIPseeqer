#!/usr/bin/perl

use lib "$ENV{HOME}/PERL_MODULES";

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();


my $n = $ta->getNbColumns();
if ($n == 0) {
  print scalar(@$a_ref) . " rows but no columns !\n";
  exit;
}
for (my $i=0; $i<$n; $i++) {
    
    my $a_ref = $ta->getColumn($i);

    print join("\t", @$a_ref); print "\n";
    

}
