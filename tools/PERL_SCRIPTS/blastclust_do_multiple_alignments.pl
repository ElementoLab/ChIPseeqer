BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;

my $ta = Table->new;
$ta->setDelim(" ");
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my $c = 0;
foreach my $r (@$a_ref) {
  
  print "Processing $c ... ";
  my $f = "TMP/c$c";
  open OUT, ">$f.txt";
  foreach my $s (@$r) {
    print OUT "$s\n";
  }
  close OUT;
  
  system("perl ~/PERL_MODULES/SCRIPTS/fasta_get_sequences_from_fasta_file.pl $f.txt all.fa > $f.fa.t");
  
  system("perl ~/PERL_MODULES/SCRIPTS/fasta_modify_names_for_clustalw.pl $f.fa.t > $f.fa");

  system("clustalw $f.fa");

  print "Done.\n";
  $c ++;
}

