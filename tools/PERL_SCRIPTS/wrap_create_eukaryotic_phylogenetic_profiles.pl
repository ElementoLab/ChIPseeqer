
BEGIN{ $home = `echo \$HOME`; chomp $home}

use lib "$home/PERL_MODULES";
use Sets;
use PBS;

my $a_ref_files = Sets::getFiles("xa*");
my $pwd = `pwd`; $pwd =~ s/\n//;

foreach my $f (@$a_ref_files) {
  
  my $pbs = PBS->new;
  
  my $prj_dir = "$home/PROJECTS/EUKARYOTIC_PHYLOGENETIC_PROFILES";
  
  $pbs->addCmd("cd $pwd");
  $pbs->addCmd("perl $home/PERL_MODULES/SCRIPTS/create_eukaryotic_phylogenetic_profiles.pl PROTEOMES/Homo_sapiens.NCBI35.sep.pep.fa PROTEOMES/list_proteomes.txt $f > $f.OUT");
  
  $pbs->submit;

}

