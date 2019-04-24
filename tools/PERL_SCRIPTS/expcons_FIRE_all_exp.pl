BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


my $root = "/home/olly/PROGRAMS/EXPCONS/ALN/yeast_gasch_IclustPos.txt.summary_EXPCON";
my @slen = (50, 125, 250, 500);

for (my $i=1; $i<26; $i++) {
  
  foreach my $l (@slen) {
    system("cp $root/$i.pro $root/$i.pro.$l");    
    system("perl fire.pl --expfile=$root/$i.pro.$l --exptype=continuous --species=yeast --fastafile_dna=$root/$i\_$l.seq --dorna=0 --dodnarna=0");    
  }
  
}
