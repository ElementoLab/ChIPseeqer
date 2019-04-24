use lib qw(/home/olly/PERL_MODULES);

use Sets;
use Table;
use Fasta;
use strict;

if (!$ARGV[1]) {
	die "Usage : prg seq kmers\n";

}


#
#  load the miRNAs

my %FPKMERS = ();
my %KMERS = ();
my %MIRNAS = ();
my $fa = Fasta->new;
$fa->setFile($ARGV[0]);
while (my $a_seq = $fa->nextSeq()) {
    my ($n, $s) = @$a_seq;

    $MIRNAS{ $n } = $s;

    $s =~ s/u/t/g;
    $s = uc($s);
    #my $c       = Sets::getComplement($s);
    my $l       = length($s);
    
    for (my $i=0; $i<=2; $i++) {

	my $k7 = Sets::getComplement(substr($s, $i, 7));
	my $k8 = Sets::getComplement(substr($s, $i, 8));
	my $k9 = Sets::getComplement(substr($s, $i, 9));
	push @{ $FPKMERS{ $k7 } }, $n;
	push @{ $FPKMERS{ $k8 } }, $n;
	push @{ $FPKMERS{ $k9 } }, $n;
	
    }
    
    for (my $i=3; $i<=$l-7; $i++) {
	my $k7 = Sets::getComplement(substr($s, $i, 7));
	push @{ $KMERS{ $k7 } }, $n;
    }

    for (my $i=3; $i<=$l-8; $i++) {
	my $k8 = Sets::getComplement(substr($s, $i, 8));
	push @{ $KMERS{ $k8 } }, $n;
    }

    for (my $i=3; $i<=$l-9; $i++) {
	my $k9 = Sets::getComplement(substr($s, $i, 9));
	push @{ $KMERS{ $k9 } }, $n;
    }


}


my $ta = Table->new;
#$ta->setDelim(" ");
$ta->loadFile($ARGV[1]);

my $a_ref = $ta->getArray();

my $only5       = $ARGV[2];
my $fulldisplay = 1;

my $i = 1;
foreach my $r (@$a_ref) {

    print "$r->[0]\t$r->[1]\n";

    if (defined($FPKMERS{ $r->[0] }) && defined($FPKMERS{ $r->[1] })) {
	print "$r->[0]\t$r->[1]\t";
	foreach my $m (@{ $FPKMERS{ $r->[0] } }) {
	    print "  $i: $r->[0]\t$m\t$MIRNAS{$m}\t0\t";
	}
	foreach my $m (@{ $FPKMERS{ $r->[1] } }) {
	    print "  $i: $r->[1]\t$m\t$MIRNAS{$m}\t0\t";
	}
	
	print "\n";
    }
    

    #foreach my $m (@{ $KMERS{ $r->[0] } }) {
#	print "$i: $r->[0]\t$m\t$MIRNAS{$m}\t3\n";
#    }

    $i++;
}
    
