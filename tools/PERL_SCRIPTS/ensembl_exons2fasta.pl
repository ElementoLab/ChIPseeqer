#!/usr/bin/perl

use lib qw(/home/olly/PERL_MODULES);

use Fasta;

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);
my $cnt = 0;
while (my $a_ref = $fa->nextSeq()) {

    my ($n, $s) = @$a_ref;

    my ($nn) = $n =~ /\|(.+?)\ /;

    if ($SEQ{$nn}) {
	$SEQ{$nn} .= "N" x 20;
	$SEQ{$nn} .= $s;
    } else {
	$SEQ{$nn}  = $s;
    }

    $cnt ++;

    #last if ($cnt == 100);
}


while (my ($k, $v) = each(%SEQ)) {

    print ">$k\n$v\n\n";
}
