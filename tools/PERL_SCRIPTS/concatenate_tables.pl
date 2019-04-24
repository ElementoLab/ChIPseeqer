#!/usr/bin/perl

use lib qw(/home/elemento/PERL_MODULES);
use Sets;
use Table;


my $ta = Table->new;
$ta->loadFile($ARGV[0]) or "Cannot load matrix1 file\n";
my $a_ref1 = $ta->getArray();

$ta->loadFile($ARGV[1]) or "Cannot load matrix2 file\n";
my $a_ref2 = $ta->getArray();

my $n1 = scalar(@$a_ref1);
my $n2 = scalar(@$a_ref2);

if ($n1 != $n2) {
    die "Problem : the two tables should have the same number of rows ..\n";
}


for (my $i=0; $i<$n1; $i++) {

    print join("\t", @{ $a_ref1->[$i] });
    
    my $m = shift @{ $a_ref2->[$i] };

    if ($m ne $a_ref1->[$i]->[0]) {
	die "Problem : the two index are different ..\n";
    }
    
    print "\t";

    print join("\t", @{ $a_ref2->[$i] });
	
    print "\n";
}
