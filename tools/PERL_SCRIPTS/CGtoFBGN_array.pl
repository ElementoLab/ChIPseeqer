use lib qw(/home/olly/PERL_MODULES);
use Fasta;
use Table;

if (scalar(@ARGV) == 0) {
	die "usage : prg fasta CGtoFB 1\n";
}

my $t = Table->new;
if ($ARGV[2]) {
 $t->setUC(1);
}

$t->loadFile($ARGV[1]);
my $h_ref = $t->getIndexColumnsKV(0, 1);

if ($ARGV[2]) {
 $t->setUC(0);
}

$t->loadFile($ARGV[0]);
my $a_ref = $t->getArray();

open LOG, ">log-CGtoFBGN_array.txt";


if ($ARGV[2]) {
  my $s = shift @$a_ref;
  print join("\t", @$s); print "\n";
}


foreach my $r (@$a_ref) {
    my $rr = $r->[0];
    #$rr =~ s/\.1$//;
    #$rr = uc($h_ref->{$rr});

    $rr = $h_ref->{$rr};

    if ($rr ne "") {
	$r->[0] = $rr;
	print join("\t", @$r); print "\n";
    } else {
	print LOG "no FBGN for \"$r->[0]\"\n";
    }
}

close LOG;

