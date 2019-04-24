use lib qw(/home/olly/PERL_MODULES);
use Fasta;


my $f = Fasta->new;

$f->setFile($ARGV[0]);

while (my $a_ref = $f->nextSeq()) {

    my ($name, $seq) = @$a_ref;
    
    #next if (length($seq) < 50);
    
    

    $name =~ s/\.1\|.+$//;
    
    if ($ARGV[1] > 0) {
	$seq = substr($seq, 0, $ARGV[1]);
    }

    print ">$name\n$seq\n\n";
    
}
