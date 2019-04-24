#!/usr/bin/perl

use lib qw(/home/olly/PERL_MODULES /home/olly/PERL_MODULES/Hypergeom /home/olly/PERL_MODULES/Hypergeom/blib/lib /home/olly/PERL_MODULES/Hypergeom/blib/arch/auto/Hypergeom);

use ScanACE;
use Table;
use GO_func;
use Yeast;
use Sets;
use Library_Motif_RE;
use GO_categories;
use DataFiles;
use strict;
use Getopt::Long;
use File::Copy;


my $m = GO_func->new;
$m->setSource("GO", "WORM");
$m->setTotalNbORFS(19873);
$m->setPvalueThreshold(0.0001);



$m->setORFset(Sets::readSet($ARGV[0]));
$m->setVerbose(0);
my $a_ref_func = $m->getFunctionalContent;
foreach my $r (@$a_ref_func) {
    print  sprintf("%3.2e\t%s\n", $r->{PVALUE}, $r->{TEXT});
    
} 


