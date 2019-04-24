use lib qw(/home/olly/PERL_MODULES);

require 'Sequence.pm';
use Sets;
use Table;
use strict;

my $s = Sequence->new;
#$s->setVerbose(1);


my $n = shift @ARGV;
foreach my $f (@ARGV) {
    $s->setBlastDB($f);
    my $seq = $s->getSequenceFromBlastDB($n, 0,0);

    if ($seq) {
	print ">$n $f\n$seq\n\n";
    }
}

