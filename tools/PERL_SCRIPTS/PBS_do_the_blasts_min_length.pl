BEGIN{ $home = `echo \$HOME`; chomp $home}

use lib "$home/PERL_MODULES";
use Sets;
use PBS;

my $a_ref_files = Sets::readSet($ARGV[0]);

my $pwd = `pwd`; $pwd =~ s/\n//;

my $cnt = 0;
foreach my $f (@$a_ref_files) {
  
  my $pbs = PBS->new;
  my $dir = $f; $dir =~ s/\/genome\.faa//g;
  $pbs->addCmd("cd $pwd");
  $pbs->addCmd("setenv PATH \${PATH}:/home/elemento/PERL_MODULES/PROGRAMS/BLAST/bin");
  $pbs->addCmd("perl $home/PERL_MODULES/SCRIPTS/do_the_blasts_min_length.pl $f $ARGV[0] > $dir/profiles.txt");
  
  $pbs->submit;
  #$pbs->print;
  $cnt ++;

  last if ($cnt == 10);
}

