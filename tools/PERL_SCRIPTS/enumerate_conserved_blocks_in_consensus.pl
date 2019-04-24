# get a set of files as input

use lib qw(/home/olly/PERL_MODULES);




use  ClustalW;
use Sets;


my $a_ref_files = Sets::readSet($ARGV[0]);

my %H = ();
foreach my $f (@$a_ref_files) {
    
    open IN, $f;
    my $l = <IN>;
    chomp $l;
    my @a = split /\-+/, $l;
    
    foreach my $r (@a) {
	if (length($r) >= 6) {

	    my $ln = length($r);
	    
	    my $maxlen = Sets::min(15, $ln);
	    for (my $len=6; $len<=$maxlen; $len++) {
		for (my $i=0; $i<length($r)-$len+1; $i++) {
		    my $rs = substr($r, $i, $len);

		    if (!defined($H{$rs}) && !defined($H{ Sets::getComplement($rs) })) { 
			$H{ $rs } = 1; #print "Add $rs (from $r)\n";
		    }
		}
	    }

	    #&& (length($r) <= 15)) {
	    #if (!defined($H{$r}) && !defined($H{ Sets::getComplement($r) })) { 
	    #	$H{ $r } = 1;
	    #}
	}
    }

    close IN;
    
    #$cnt ++;

    #last if ($cnt == 10);
}


foreach my $r (keys(%H)) {
    print "$r\n";
}
