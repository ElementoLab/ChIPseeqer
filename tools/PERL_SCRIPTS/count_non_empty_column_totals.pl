#!/usr/bin/perl
#
#  usage : prg thr 7mers 8mers 9mers ..
#
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use Table;
use strict;

my $ta = Table->new;

$ta->loadFile($ARGV[0]);

my $a_ref = $ta->getArray();
my @CNT = ();
foreach my $r (@$a_ref) {
    
    my $n = scalar(@$r);
    
    for (my $i=0; $i<$n; $i++) {
	my $s = $r->[$i]; $s =~ s/\ //g;
	$CNT[$i] ++ if ($s ne "");
    }

}


print join("\t", @CNT) . "\n";
