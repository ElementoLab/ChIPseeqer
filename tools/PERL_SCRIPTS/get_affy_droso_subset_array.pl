use lib qw(/home/olly/PERL_MODULES);


use Table;
use Sets;
use Getopt::Long;
use Drosophila;
#use strict;

my $type = "u";
my $pma     = undef;
my $outdata = undef;
my $logt    = 0;
GetOptions ('conds=s'     => \$conds,
	    'file=s'      => \$file,
            'genesonly=s' => \$genesonly,
	    'chr=s'       => \$chr,
	    'logt=s'      => \$logt,	    
	    'type=s'      => \$type,
	    'pma=s'       => \$pma,    # minimal number of P
	    'outdata=s'   => \$outdata,
	    't=s'         => \$t
	    );

my @a_conds = split /\,/, $conds;

my $dr = undef;
if (defined($chr)) {
    $dr = Drosophila->new();
}

my $ta = Table->new;

$ta->loadFile($file);


my $a_ref = $ta->getArray();

my @G = ();
my $r = shift @$a_ref;
my @firstline = ("");

foreach my $c (@a_conds) {
    push @firstline, $r->[$c];
}
    my $l = join("\t", @firstline);
$l =~ s/\ //g;
$l =~ s/\#//g;
print "$l\n";


my @COUNTS = ();
foreach my $r (@$a_ref) {

    if (defined($chr)) {
	if (!$dr->isOnChromosome($r->[0], $chr)) {
	    next;
	}
    }

    #print join("\t", @$r) . "\n";

    #print "TTTTT\n";
    my $a_ref1 = [];
    my $a_ref1_PMA = [];
    foreach my $c (@a_conds) {
	#print "$r->[$c]\n";
	push @$a_ref1, $r->[$c];
	push @$a_ref1_PMA, $r->[$c+1];

	push @{ $COUNTS[ $c ] }, $r->[$c];
    }
    
    #
    #   present genes
    #

    my $pre1 = 0;
    
    if (defined($pma)) {

	foreach my $p (@$a_ref1_PMA) {
	    $pre1 ++ if ($p eq "P");
	}

	if ($pre1 < $pma) {
	    next;
	}
    }



 
    print "$r->[0]\t";
    if ($logt == 1) {
	@$a_ref1 = map(log, @$a_ref1);
	    
    }
    print join("\t", @$a_ref1); print "\n";
    
    


    
}

