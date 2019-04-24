#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use PBS;
my $pwd  = `pwd`; $pwd =~ s/\n//;

for (my $i=1; $i<=100; $i++) {
  

  my $pbs = PBS->new;
  $pbs->setWallTime("24:00:00");
  $pbs->setScriptName("script.run.AlignaCE.$i");
  $pbs->addCmd("cd $pwd");

  my $ff = "FIRE_DATA/YEAST/EXPFILES/yeast_gasch_IclustPos.txt.$i\_FIRE/DNA/yeast_gasch_IclustPos.txt.$i";
  print "Analyzing $ff\n";
  my $todo = "perl clusters_run_alignace.pl --clusters=$ff --fastafile=FIRE_DATA/YEAST/SEQUENCES/yeast_u_600_0.fa > $ff.AlignACE3";

  $pbs->addCmd($todo);
  #$pbs->execute();
  $pbs->submit();

}

