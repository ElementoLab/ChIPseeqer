use lib qw(/home/olly/PERL_MODULES);

use lalign;
use Fasta;
use Sets;


my $fa = Fasta->new;
$fa->setFile($ARGV[0]);
my $a_ref = $fa->nextSeq();
my ($n1, $s1) = @$a_ref;
$fa->dispose;

my $fa = Fasta->new;
$fa->setFile($ARGV[1]);

my $outfile  = $ARGV[2];

my @POS = ();
my $i = 1;
while (my $a_ref = $fa->nextSeq()) {
    my ($n, $s) = @$a_ref;
    
    # exec lalign
    my $file = Sets::getTempFile("/tmp/toto");
    my $outfile = Sets::getTempFile("/tmp/out");

    #print "$file\n";
    $fa->writeSeq($file, $n, $s);

    #my $todo = "/home/olly/MPLPC18/DOCTORAT/PROGRAMMES/modified/fasta2/lalign -n -f -30 $ARGV[0] $file > $outfile 2> /dev/null";
    #print "$todo\n";
    #system($todo);

    my $la = lalign->new;
    $la->setFile($outfile);
    $la->setMaxEvalue(0.1);
    #$la->process;
    
    
    $la->processUsingBlast($ARGV[0], $file);


    my $a_ref_res = $la->getResults;

    foreach my $h (@$a_ref_res) {
	
	my $a_ref_pos = $h->{CONSERVED_POS};

	foreach my $p (@$a_ref_pos) {
	    push @{ $POS[ $p ] }, $i if !Sets::in_array($i, @{ $POS[ $p ] }); # ++;
	}
	
	
    }

    $i++;
    
    #last if ($i == 3);


    unlink $file;
    unlink $outfile;
}


my $cnt = $i-1;

my @a = split //, $s1;

my $w = 100;
for (my $i=0; $i<length($s1); $i += $w) {

    my $p = $i+1;
    
    for (my $j=$cnt; $j>=0; $j--) {
	for (my $k=$i; $k<$i+$w; $k++) {
	    
	    if (defined($POS[$k]->[$j])) {
		print $POS[$k]->[$j];
	    } else {
		print " ";
	    }
	}
	print "\n";
    }

    print substr($s1, $i, $w); print "\n\n";

    #print $p . "\t" . $a[$i]; print " ";
    
    #print join("", @{ $POS[$i] }); print "\n";
}


if (defined($outfile)) {

    open OUT, ">$outfile";
    for (my $i=0; $i<length($s1); $i += 1) {
	print OUT $a[$i]; print OUT "\t";
	if ($POS[ $i ] > 0) {
	    print OUT join("", @{ $POS[$i] }); print OUT "\n";
	} else {
	    print OUT "\n";
	}
    }
    close OUT;
    
   
}

    #if ($POS[ $i ] > 0) {
#	print $POS[ $i ];
#    } else {
#	print " ";
#    }
#    print "\n";
    #print "\n" if (($i % 50) == 0);
#}
#print "\n";

#for (my $i=0; $i<length($s1); $i++) {
#    print $a[$i];
    #print "\n" if (($i % 50) == 0);
#}
#print "\n";
