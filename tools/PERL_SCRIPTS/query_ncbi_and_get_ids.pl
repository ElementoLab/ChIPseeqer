
use XML::DOM;
use LWP::Simple;
use strict;

open LOG, ">>log.txt";

my $q = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucleotide&term=%22$ARGV[0]%22&email=elemento\@princeton.edu";
my $r      = get($q);
print LOG "$q\n$r\n";
close LOG;


my $parser = new XML::DOM::Parser;
my $doc    = $parser->parse($r);


my $nodes  = $doc->getElementsByTagName("Id");
my $n      = $nodes->getLength;

#print "$n\n";

print "$ARGV[0]";
for (my $i = 0; $i < $n; $i++) {
    my $id = $nodes->item ($i)->getChildNodes->item(0)->getData;
    print " $id";
}

print "\n";

$doc->dispose;

sleep(1);


