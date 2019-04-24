
use XML::DOM;
use LWP::Simple;
use strict;



my $parser = new XML::DOM::Parser;
my $doc    = $parser->parsefile($ARGV[0]);


my $protein_node     = $doc->getElementsByTagName("protein")->item(0);


my $id               = $protein_node->getAttributeNode ("id")->getValue;


my $title_node       = $protein_node->getElementsByTagName("title")->item(0); 


my $title            = $title_node->getChildNodes->item(0)->getData;


print "$id\t$title\n";

$doc->dispose;



