BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;

open IN, $ARGV[0];
while (my $l = <IN>) {
  chomp $l;
  
  next if ($l eq "");

  

  if ($l =~ /^\>/) {

    print "$l\n";
    
    my $lh = <IN>; $lh = uc($lh);
    my $lm = <IN>; $lm = uc($lm);
    my $lr = <IN>; $lr = uc($lr);
    my $lc = <IN>; $lc = uc($lc);
    
    #$lh   =~ s/\#/N/g;
    $lh   =~ s/\-//g;

    print "$lh\n";
    
    #<STDIN>;
  }
  
 

}
close IN;
