
use XML::DOM;
use LWP::Simple;
use strict;



my $parser = new XML::DOM::Parser;
my $doc    = $parser->parsefile($ARGV[0]);


my $id     = $doc->getElementsByTagName("protein")->item(0)->getAttributeNode ("id")->getValue;
my $cdna   = $doc->getElementsByTagName("cdna")->item(0)->getChildNodes->item(0)->getData;

$cdna =~ s/ //g; 
$cdna =~ s/\n//g; 


print ">$id\n$cdna\n\n";

$doc->dispose;



