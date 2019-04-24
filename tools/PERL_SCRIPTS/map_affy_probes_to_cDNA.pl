BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Fasta;
use Sets;
use MyBlast;
use strict;



my $fa = Fasta->new;

my $mb = MyBlast->new;
$mb->setBlastProgram("blastn");

$mb->setDatabaseDatabase($ARGV[1]);
$mb->setNbProcessors(2);
$mb->setEvalueThreshold("1e-5");
$mb->setQueryStrand(1);
$mb->setVerbose(0);
    
    
#
# go thru all the sequences, align them one by one
#

my $file = Sets::getTempFile("toto");

open IN, $ARGV[0];
my $l = <IN>;
while (my $l = <IN>) { 
    chomp $l;
    
    my @a = split /\t/, $l;

    my $n = $a[0];
    my $s = $a[4];
    
    $fa->writeSeq($file, $n, $s);
    $mb->setQueryFile($file);
    
    my $a_ref = $mb->blastallMultiple;

    
    #print "got " . scalar(@$a_ref) . " hits\n";

    foreach my $hit (@$a_ref) { 
	
	my $hsps = $hit->{"HSPS"};
	
	foreach my $r (@$hsps) {
	    
	    my $s1 = $r->{"ALIGNLEN"};
	    my $e1 = $r->{"IDENTITY"};
	    my $ev = $r->{"EVALUE"};
	    my $st = $r->{"DFRAME"};
	    
	    
	    #$n =~ s/\#/\t/;
	    
	    next if ($st == -1);

	    print "$n\t$hit->{HIT_NAME}\t$e1\t$ev\t$st\n";
	}
    }

    #<STDIN>;
}

unlink $file;
