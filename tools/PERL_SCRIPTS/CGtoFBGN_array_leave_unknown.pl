use lib qw(/home/olly/PERL_MODULES);
use Fasta;
use Table;

if (scalar(@ARGV) == 0) {
	die "usage : prg array CGtoFB\n";
}

my $t = Table->new;
$t->setUC(1);

$t->loadFile($ARGV[1]);
my $h_ref = $t->getIndexColumnsKV(0, 1);



$t->loadFile($ARGV[0]);
my $a_ref = $t->getArray();

foreach my $r (@$a_ref) {
    $r->[0] =~ s/\.1$//;
    if (defined($h_ref->{$r->[0]})) {
	$r->[0] = uc($h_ref->{$r->[0]});
    }

    #if ($r->[0] ne "") {
    print join("\t", @$r); print "\n";
    #} 
}

