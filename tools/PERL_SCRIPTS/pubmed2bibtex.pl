#!/usr/bin/perl

use lib qw(/home/olly/PERL_MODULES);
use MyBiblio;


my $b = MyBiblio->new;

$b->getArticleByPMID($ARGV[0]);

print $b->toBibTex;

