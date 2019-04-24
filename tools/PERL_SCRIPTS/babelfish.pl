#!/usr/bin/perl

use HTML::Form;
use LWP::Simple;
use LWP;

$HTML = get("http://world.altavista.com/");

@forms = HTML::Form->parse($HTML, "http://world.altavista.com/");

#print scalar(@forms) . "\n";

$f = shift @forms;
    
#@inputs = $f->inputs;

#foreach $i (@inputs) {
#    print $i->type . "\t" .  $i->name . "\t" . $i->value . "\n";
    
#}


$f->value("urltext", $ARGV[0]);

$f->value("lp", ($ARGV[1]?$ARGV[1]:"fr_en"));

my $ua = LWP::UserAgent->new;

$response = $ua->request($f->click);
    
#print $response->content;

$RES = $response->content;

$RES =~ /\<td bgcolor=white class=s\>\<div style=padding\:10px\;\>(.+?)\<\/div\>\<\/td\>/;

print "$1\n";;


    

