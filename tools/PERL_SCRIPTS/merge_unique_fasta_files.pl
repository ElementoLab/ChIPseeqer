use lib qw(/home/olly/PERL_MODULES);

use Fasta;

my %H = ();
foreach my $f (@ARGV) {
    my $fa1 = Fasta->new;
    $fa1->setFile($f);
    
    while (my $a_ref1 = $fa1->nextSeq()) {
    
	my ($n, $s) = @$a_ref1;
	
	my ($n1, $n2) = split /\s/, $n;
	
	if (!defined($H{$n1})) {
	    
	    

	    print ">$n\n$s\n\n";
	    $H{$n1} = $n2;

	} 

	  
    }
    
}



