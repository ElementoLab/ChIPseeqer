
use XML::DOM;
use LWP::Simple;
use strict;

#print "$ARGV[0]\n";

my $parser = new XML::DOM::Parser;
my $doc    = $parser->parsefile($ARGV[0]);


my $id     = $doc->getElementsByTagName("protein")->item(0)->getAttributeNode ("id")->getValue;

my $interactionListNodes   = $doc->getElementsByTagName("interactionList");
my $ni                     = $interactionListNodes -> getLength;

exit if ($ni == 0);

my $interactionList   = $interactionListNodes->item(0);



my $interactions      = $interactionList->getElementsByTagName("interaction");

my $n                 = $interactions -> getLength;
#print $n;

for (my $i = 0; $i < $n; $i++) {
    

    my @PMID = ();
    my $bibref      = $interactions->item($i)->getElementsByTagName("bibref")->item(0);
    
    next if !defined($bibref);

    foreach my $priref ($bibref->getElementsByTagName("primaryRef")) {
	my $db = $priref->getAttributeNode ("db")->getValue; 
	my $id = $priref->getAttributeNode ("id")->getValue; 

	if ($db eq "PubMed") {
	    push @PMID, $id;
	}

	

    }
    
    my @INT = ();
    foreach my $p  ($interactions->item($i)->getElementsByTagName("proteinInteractorRef")) {

	my $ref = $p->getAttributeNode ("ref")->getValue; 
	$ref =~ s/ID_//;
	
	push @INT, $ref;
    }
    
    my $np = scalar(@INT);
    
    for (my $j=0; $j<$np-1; $j++) {
	for (my $k=$j+1; $k<$np; $k++) {
	    print "$INT[$j]\t$INT[$k]";
	    print "\t"; print join("/", @PMID); print "\n";
	}
    }

    

}


$doc->dispose;



