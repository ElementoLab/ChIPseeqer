#!/usr/bin/perl

use HTML::Form;
use LWP::Simple;
use LWP;

$HTML = get("http://db.yeastgenome.org/cgi-bin/SGD/locus.pl?locus=" . $ARGV[0]);


#$HTML =~ /\<tr\>\<td valign=\"top\" nowrap=\"\" bgcolor=\"\#CCCCFF\"\>Description \<\/td\>\<td valign=\"top\"\>(.+?)\<\/td\>\<\/tr\>/;

$HTML =~ /\<tr\>\<td valign=\"top\" width=\"160\" nowrap=\"\" bgcolor=\"\#CCCCFF\"\>Description\<\/td\>\<td width=\"3\"\> \<\/td\>\<td valign=\"top\"\>(.+?)\<\/td\>\<\/tr\>/;


print "$ARGV[0]: $1\n";
    

    

