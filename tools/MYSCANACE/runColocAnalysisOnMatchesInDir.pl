#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";

use Sets;
use PBS;

my $a_ref_f = Sets::getFiles("$ARGV[0]/*.matches.txt");

foreach my $f1 (@$a_ref_f) {
  
  foreach my $f2 (@$a_ref_f) {
    
    next if ($f1 eq $f2);
   
    my $pbs = PBS->new;
    
    my $ff1 = Sets::filename($f1);
    my $ff2 = Sets::filename($f2);

    $pbs->setScriptName("script.pbs.$ff1.$ff2");
    $pbs->setWallTime("24:00:00");
    my $todo = "$ENV{MYSCANACEDIR}/EvaluateMotifColocalization $ENV{HOME}/PROGRAMS/PLASMODIUM/GENOMES/sequences/6.0/pf_u_2000_0_masked.fa $f1 $f2";
    $pbs->addCmd($todo);

    $pbs->submit;
    #system($todo);
    
    
  }
	

}
