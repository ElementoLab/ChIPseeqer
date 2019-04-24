use lib qw(/home/elemento/PERL_MODULES);
use Table;
use Sets;

use Fasta;
use strict;




my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getColumn(0);

#my $a_ref    = Sets::readSet($ARGV[0]);

my $seqfile  = $ARGV[1];


my %H = ();
foreach my $r (@$a_ref) {
    $H{ $r } = [];
    
    my $k = length($r);
    die "k!=7,8, please modify the program ($r) \n" if (length($r) > 9);
}



my $fa = Fasta->new;
$fa->setFile($ARGV[1]);

while (my $a_ref = $fa->nextSeq()) {
    my ($n, $s) = @$a_ref; 
    
    my $l = length($s);
    for (my $i=0; $i<$l-7; $i++) {
	
	my $ss = substr($s, $i, 7);
	if (defined($ss) && ($ss !~ /\-/)) {
	    push @{ $H{ $ss } }, $n if (!Sets::in_array($n, @{ $H{ $ss } })); 
	}

	my $ss = substr($s, $i, 9);
	if (defined($ss) && ($ss !~ /\-/)) {
	    push @{ $H{ $ss } }, $n if (!Sets::in_array($n, @{ $H{ $ss } })); 
	}

	
	
    }
}



foreach my $r (@$a_ref) {
  
  if (scalar(@{ $H{ $r } }) > 0) {
    
    open OUT, ">$ARGV[2]/$r.txt";
    print OUT join("\n", @{ $H{ $r } }); print OUT "\n";
    close OUT;
  }
}
