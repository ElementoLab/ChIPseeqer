#!/usr/bin/perl 

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use Motif;


my $mwm = Motif->new;
$mwm->readScanACEMotif($ARGV[0]);


print $mwm->getSimpleMatrix();

