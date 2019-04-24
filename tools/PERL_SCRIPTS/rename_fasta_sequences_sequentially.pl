use lib qw(/home/olly/PERL_MODULES);
use Fasta;


my $f = Fasta->new;

$f->setFile($ARGV[0]);

my $cnt = 1;
while (my $a_ref = $f->nextSeq()) {

    my ($name, $seq) = @$a_ref;
    
    #next if (length($seq) < 50);
    
   $name =~ s/\ .+$//; 

    $name .= ".$cnt";
    

    print ">$name\n$seq\n\n";
    

    $cnt++;
}
