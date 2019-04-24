# perl ~/PERL_MODULES/SCRIPTS/expression_get_U95A_probe_NM_annotation.pl GPL91.txt  
# hg17_up2k_dwn2k_aligned.nmbased_humanonly.NMids.txt > GPL91_probe_NM_mapping.txt

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;
use strict;

my $h_ref = Sets::getIndex($ARGV[1]);

open IN, $ARGV[0];
while (my $l = <IN>) {
  
  next if ($l =~ /^[\^\#\!]/); 
    
  my @a  = split /\t/, $l, -1;
  
  next if ($a[10] eq "");

  $a[10] =~ s/ \/\/\/ /\t/g;
  
  my @b = split /\t/, $a[10];

  my @c = ();
  foreach my $bb (@b) {
    if (defined($h_ref->{$bb})) {
      push @c, $bb;
    }
  }

  if (@c != 0) {
    print "$a[0]\t" . join("\t", @c) . "\n";
  }
  
  

}
