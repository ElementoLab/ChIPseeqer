use lib qw(/home/olly/PERL_MODULES);

require 'Sequence.pm';


if (scalar(@ARGV) < 3) {
    die "Usage : getseqsfromlist.pl DB1 DB2 DB3 ORT MIN MAX\n";
}


my $thedbfile1  =  $ARGV[0];
my $thedbfile2  =  $ARGV[1];
my $thedbfile3  =  $ARGV[2];
my $thegenefile =  $ARGV[3];
my $min         = ($ARGV[4]?$ARGV[4]:100);
my $max         = ($ARGV[5]?$ARGV[5]:1000);


open OUT1, ">upstream_sp1.fasta";
open OUT2, ">upstream_sp2.fasta";
open OUT2, ">upstream_sp3.fasta";


open INCH, $thegenefile;
while (my $l = <INCH>) {
    
    chomp $l;

    my @a = split /\t/, $l;
	    
    my $o_seq = Sequence->new;

#$o_seq->setVerbose(1);
    $o_seq->setBlastPath("/home/olly/COMPARATIVE_YEAST/ORFS/BLAST_REPORTS/BINDINGMAP/PROGRAMS/BLAST");

    $o_seq->setBlastDB($thedbfile1);    
    my $s1 = $o_seq->getSequenceFromBlastDB($a[0], 0, 0);
    
    $o_seq->setBlastDB($thedbfile2);    
    my $s2 = $o_seq->getSequenceFromBlastDB($a[1], 0, 0);

    $o_seq->setBlastDB($thedbfile3);    
    my $s3 = $o_seq->getSequenceFromBlastDB($a[2], 0, 0);

    if ((length($s1) > $min) && (length($s2) > $min) && (length($s3) > $min)) {

	if ($max) {
	    $s1 = substr($s1, (length($s1)<$max?0:length($s1) - $max), $max);
	    $s2 = substr($s2, (length($s2)<$max?0:length($s2) - $max), $max);
	    $s3 = substr($s3, (length($s3)<$max?0:length($s3) - $max), $max);

	}

	print OUT1 ">$a[0]\n$s1\n\n";
	print OUT2 ">$a[1]\n$s2\n\n";
	print OUT3 ">$a[2]\n$s3\n\n";

    }	
    
}
    
    
close INCH;
