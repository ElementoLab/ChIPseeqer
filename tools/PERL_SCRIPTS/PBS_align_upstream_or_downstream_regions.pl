my $home = undef;

BEGIN{ $home = `echo \$HOME`; chomp $home; }

use lib "$home/PERL_MODULES";
use Sets;
use PBS;
use strict;

system("mkdir $ARGV[1]") if (! -e "$ARGV[1]");
my $a_ref_genes = Sets::readSet($ARGV[0]);
my $pwd = `pwd`; $pwd =~ s/\n//;

my $cnt = 0;
foreach my $f (@$a_ref_genes) {
  
  my $pbs = PBS->new;
    
  $pbs->addCmd("cd $pwd");
  $pbs->addCmd("setenv DIALIGN2_DIR /home/elemento/PERL_MODULES/PROGRAMS/DIALIGN");
  $pbs->addCmd("perl $home/PERL_MODULES/SCRIPTS/align_upstream_or_downstream_regions.pl $f U 2000 DATA/ortholog_table.txt");
  $pbs->addCmd("mv $f.ali $ARGV[1]/");  
  $pbs->addCmd("rm $f.seq");  
  $pbs->submit;
  
  $cnt ++;

  #last if ($cnt == 10);
}

