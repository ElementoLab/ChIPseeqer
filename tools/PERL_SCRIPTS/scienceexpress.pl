#!/usr/bin/perl

use HTML::Form;
use LWP::Simple;
use LWP;

my $HTML = get("http://www.sciencemag.org/sciencexpress/recent.shtml");



while ($HTML =~ /\<STRONG\>(.+?)\<\/STRONG\>\<DD\>/ig) {
    my $l = $1;
    $l =~ s/<.+?>//g;
    print "$l\n";
}

    

