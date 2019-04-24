use lib qw(/home/olly/PERL_MODULES);
use Fasta;
use Table;

if (scalar(@ARGV) == 0) {
	die "usage : prg fasta CGtoFB\n";
}

my $t = Table->new;
$t->setUC(1);

$t->loadFile($ARGV[1]);
my $h_ref = $t->getIndexColumnsKV(0, 1);


my $f = Fasta->new;

$f->setFile($ARGV[0]);

while (my $a_ref = $f->nextSeq()) {
    my ($name, $seq) = @$a_ref;
    $name =~ s/\.1$//;
    $name = uc($h_ref->{$name});
    
    if ($name ne "") {
	print ">$name\n$seq\n\n";
    } 
}
