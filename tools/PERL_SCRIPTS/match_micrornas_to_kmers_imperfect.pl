use lib qw(/home/olly/PERL_MODULES);

use Sets;
use Table;
use Fasta;
use strict;

if (!$ARGV[1]) {
	die "Usage : prg seq kmers\n";

}

my $ta = Table->new;
$ta->loadFile($ARGV[1]);

my $a_ref = $ta->getArray();

my $only5       = $ARGV[2];
my $fulldisplay = 1;

my $i = 1;
foreach my $r (@$a_ref) {
    
    
    my $fa = Fasta->new;
    $fa->setFile($ARGV[0]);

    my $match = 0;
    my $rank  = 1;
    while (my $a_seq = $fa->nextSeq()) {
	
	my ($n, $s) = @$a_seq;

	my $ss = substr($s, 0, 8);

	$ss =~ s/t/u/g;
	
	$ss = uc($ss);
	$ss =~ s/U/T/g;

	$ss = Sets::getComplement($ss);

	#for (my $j=0; $j<=0; $j++) { #  uhh ?
	my $j = 0;

	my $kmer = substr($r->[0], 0, length($r->[0])-$j);
	my $km   = $kmer;
	
	$km =~ s/N/\./g;

	#
	#  but here we are interested only in 
	#
	my $len = length($km);
	for (my $j=0; $j<$len; $j++) {
	    my @a = split //, $km;
	    $a[$j] = "[^$a[$j]]";
	    my $km_mod = join("", @a); #print "rem=$km_mod($km)\n";
	    my $p = undef;
	    if ($ss =~ /($km_mod)/) {
		my $pat = Sets::getComplement($1); $pat = lc($pat); $pat =~ s/t/u/g; $pat = reverse($pat);

		my $i1 = substr($r->[0], $j, 1);  $i1 = lc($i1); $i1 =~ s/t/u/g; 
		my $i2 = substr($pat   , $j, 1); 

		print "$i: $r->[0]\t$n\t$s\t$km_mod\t$r->[0]-$pat\t($i1-$i2)\n" if (("($i1-$i2)" eq "(u-g)") || ("($i1-$i2)" eq "(g-u)")); 
	    }
	}
	
    }

    $fa->dispose();
    
    $i++;
}
