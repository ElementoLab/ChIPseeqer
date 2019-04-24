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
GetOptions ('cols=s'     => \$cols,
            'repl=s'     => \$repl,
            'file=s'      => \$file,
            'genesonly=s' => \$genesonly,
	    'chr=s'       => \$chr,
	    'logt=s'       => \$logt,
	    
	    'type=s'      => \$type,
	    'pma=s'       => \$pma,    # minimal number of P
	    'outdata=s'   => \$outdata,
	    't=s'         => \$t);

my @a_cols = split /\,/, $cols;

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
if ($outdata) {
    foreach my $c (@a_cond1) {
	push @firstline, $r->[$c];
    }
    foreach my $c (@a_cond2) {
	push @firstline, $r->[$c];
    }
    my $l = join("\t", @firstline);
    $l =~ s/\ //g;
    $l =~ s/\#//g;
    print "$l\n";
}

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
    foreach my $c (@a_cols) {
	#print "$r->[$c]\n";
	push @$a_ref1, $r->[$c];
	push @$a_ref1_PMA, $r->[$c+1];

	push @{ $COUNTS[ $c ] }, $r->[$c];
    }

    my $a_ref2 = [];
    my $a_ref2_PMA = [];
    foreach my $c (@a_cond2) {
	#print "$r->[$c]\n";
	push @$a_ref2, $r->[$c];
	push @$a_ref2_PMA, $r->[$c+1];

	push @{ $COUNTS[ $c ] }, $r->[$c];
    }
    
    #
    #   present genes
    #
    my $pre1 = 0;
    my $pre2 = 0;
    if (defined($pma)) {

	foreach my $p (@$a_ref1_PMA) {
	    $pre1 ++ if ($p eq "P");
	}

	foreach my $p (@$a_ref2_PMA) {
	    $pre2 ++ if ($p eq "P");
	}
	
	#print "$pre1\t$pre2\n";

	#if (($pre1 < $pma) && ($pre2 < $pma)) {
	if ($pre1+$pre2 < $pma) {
	    next;
	}
    }



    if ($outdata) {
	print "$r->[0]\t";
	if ($logt == 1) {
	    @$a_ref1 = map(log, @$a_ref1);
	    @$a_ref2 = map(log, @$a_ref2);
	}
	print join("\t", @$a_ref1); print "\t";
	print join("\t", @$a_ref2); print "\n";
	
    }

    
    
    #print join("\t", @$a_ref2); print "\n";

    #my $go = (Sets::average($a_ref2) - Sets::average($a_ref1)) / 
    #	(Sets::stddev($a_ref2) + Sets::stddev($a_ref1));

    my $go = undef;

    $type = lc($type);
    
    if (!$outdata) {
    
    if ($type eq "d") {
	$go = sprintf("%3.2f", Sets::average($a_ref1) / Sets::average($a_ref2));
    } else {
	$go = sprintf("%3.2f", Sets::average($a_ref2) / Sets::average($a_ref1));
    }
    
    }

    if ($go >= $t) {
	
	#print Sets::average($a_ref1);
	#print "\t";
	#print Sets::average($a_ref2);
	#print "\n";

	my @a = ($r->[0], $go);

	push @G, \@a;
	

    }
}


my $n = scalar(@COUNTS);
for (my $i=0; $i<$n; $i++) {
    #print join("\n", @{ $COUNTS[$i] });
    print Sets::median( $COUNTS[$i] ); print "\n";
}

#exit if ($outdata);

#@G = sort {$b->[1] <=> $a->[1] } @G;

#foreach my $r (@G) {
#    print "$r->[0]";
#    print "\t$r->[1]" if ($genesonly != 1);
#    print "\n";
#
#};
