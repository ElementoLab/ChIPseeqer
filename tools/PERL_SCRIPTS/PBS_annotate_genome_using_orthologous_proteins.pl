BEGIN{ $home = `echo \$HOME`; chomp $home}

use lib "$home/PERL_MODULES";
use Sets;
use PBS;

system("mkdir TMP") if (! -e "TMP");
system("split -l $ARGV[3] $ARGV[2] TMP/");
my $a_ref_files = Sets::getFiles("TMP/*");
my $pwd = `pwd`; $pwd =~ s/\n//;

foreach my $f (@$a_ref_files) {
  
  my $pbs = PBS->new;
    
  $pbs->addCmd("cd $pwd");
  $pbs->addCmd("setenv WISECONFIGDIR /home/elemento/PERL_MODULES/PROGRAMS/GENEWISE/wise2.2.0/wisecfg/");
  $pbs->addCmd("perl $home/PERL_MODULES/SCRIPTS/annotate_genome_using_orthologous_proteins.pl $ARGV[0] $ARGV[1] $f > $f.OUT");
  
  $pbs->submit;

}

