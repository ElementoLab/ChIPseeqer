use lib qw(/home/olly/PERL_MODULES);
use Sets;
use ClustalW;


my $cl = ClustalW->new;

my $a_ref = Sets::getFiles($ARGV[0]);

my %OUT = ();

foreach my $r (@$a_ref) {
    
    my $c = $r;
    $c =~ s/\.seq//;
    
    next if (! -e "$c.aln");

    

    #next if ($c ne "CG10138");
    print "$c\n";
    #my $todo = "clustalw $r -OUTORDER=INPUT";

    #print $todo;

    #system($todo) == 0 or next;

    $cl->setFile("$c.aln");

    my $h_ref_seqs = $cl->getSeqs;
    
    my $cnt = scalar(keys(%$h_ref_seqs));
    next if ($cnt != 7);


    foreach my $s (keys(%$h_ref_seqs)) {
	
	#print "$s\n";

	#die "died $c" if (length($s) == 0);
	#print "adding to $s sequence $h_ref_seqs->{$s}\n";
	$OUT{$s} .= $h_ref_seqs->{$s};

	$cnt ++;
    }

    
    #foreach my $s (keys(%OUT)) {
#	print "$s\t" . length($OUT{$s}) . "\n";
#    }

 #   print "\n\n\n\n";

}


open OUT, ">align";
foreach my $s (keys(%OUT)) {
    print OUT "$s\t$OUT{$s}\n";
}
close OUT;
