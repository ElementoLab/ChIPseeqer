#!/usr/bin/perl

use HTML::Form;
use LWP::Simple;
use LWP;

$HTML = get("http://caltech.wormbase.org/db/gene/gene?name=" . $ARGV[0]);


#$HTML =~ /\<tr\>\<td valign=\"top\" nowrap=\"\" bgcolor=\"\#CCCCFF\"\>Description \<\/td\>\<td valign=\"top\"\>(.+?)\<\/td\>\<\/tr\>/;

my $pat = "Concise Description\: \<\/th\> \<td class\=\"databody\"\>(.+?)\<\/td\>";

$HTML =~ /$pat/;


print "$ARGV[0]: $1\n";
    

    

