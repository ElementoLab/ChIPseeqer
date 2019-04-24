use lib qw(/home/olly/PERL_MODULES);


use Table;
use Sets;
use Getopt::Long;

my $type = "u";
my $pma     = undef;
my $outdata = undef;

GetOptions ('cond1=s'     => \$cond1,
            'cond2=s'     => \$cond2,
            'file=s'      => \$file,
            'genesonly=s' => \$genesonly,
	    'type=s'      => \$type,
	    'pma1=s'       => \$pma1,    # minimal number of P
	    'pma2=s'       => \$pma2,    # minimal number of P

	    'outdata=s'   => \$outdata,
	    't=s'         => \$t);

my @a_cond1 = split /\,/, $cond1;
my @a_cond2 = split /\,/, $cond2;


my $ta = Table->new;

$ta->loadFile($file);


my $a_ref = $ta->getArray();

my @G = ();
my $r = shift @$a_ref;
my @firstline = ("");
if ($outdata) {
    foreach my $c (@a_cond1) {
	push @firstline, $r->[$c];
    }
    foreach my $c (@a_cond2) {
	push @firstline, $r->[$c];
    }
    print join("\t", @firstline);
    print "\n";
}

foreach my $r (@$a_ref) {

    #print join("\t", @$r) . "\n";

    #print "TTTTT\n";
    my $a_ref1 = [];
    my $a_ref1_PMA = [];
    foreach my $c (@a_cond1) {
	#print "$r->[$c]\n";
	push @$a_ref1, $r->[$c];
	push @$a_ref1_PMA, $r->[$c+1];
    }

    my $a_ref2 = [];
    my $a_ref2_PMA = [];
    foreach my $c (@a_cond2) {
	#print "$r->[$c]\n";
	push @$a_ref2, $r->[$c];
	push @$a_ref2_PMA, $r->[$c+1];
    }
    
    #
    #   present genes
    #
    my $pre1 = 0;
    my $pre2 = 0;
    if (defined($pma1)) {

	foreach my $p (@$a_ref1_PMA) {
	    $pre1 ++ if ($p eq "P");
	}

	foreach my $p (@$a_ref2_PMA) {
	    $pre2 ++ if ($p eq "P");
	}
	
	#print "$pre1\t$pre2\n";

	#if (($pre1 < $pma) && ($pre2 < $pma)) {
	if (!(($pre1 >= $pma1) && ($pre2 >= $pma2))) {
	    next;
	}
    }



    if ($outdata) {
	print "$r->[0]\t";
	print join("\t", @$a_ref1); print "\t";
	print join("\t", @$a_ref2); print "\n";
	
    }

    
    
    #print join("\t", @$a_ref2); print "\n";

    #my $go = (Sets::average($a_ref2) - Sets::average($a_ref1)) / 
    #	(Sets::stddev($a_ref2) + Sets::stddev($a_ref1));

    my $go = undef;

    $type = lc($type);

    if ($type eq "d") {
	$go = sprintf("%3.2f", Sets::average($a_ref1) / Sets::average($a_ref2));
    } else {
	$go = sprintf("%3.2f", Sets::average($a_ref2) / Sets::average($a_ref1));
    }
    
    $t = -1000000000.0;

    if ($go >= $t) {
	
	#print Sets::average($a_ref1);
	#print "\t";
	#print Sets::average($a_ref2);
	#print "\n";

	my @a = ($r->[0], $go);

	push @G, \@a;
	

    }
}

exit if ($outdata);

@G = sort {$b->[1] <=> $a->[1] } @G;

foreach my $r (@G) {
    print "$r->[0]";
    print "\t$r->[1]" if ($genesonly != 1);
    print "\n";

};
