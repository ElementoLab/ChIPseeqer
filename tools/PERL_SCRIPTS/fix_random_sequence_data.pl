#!/usr/bin/perl
use lib qw(/home/olly/PERL_MODULES);
use Fasta;

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

while (my $a_ref = $fa->nextSeq) {
    my ($n, $s) = @$a_ref;

    if ($s eq "") {
	$s = 'N' x 2000;
    }

    print ">$n\n$s\n\n";
    
}
