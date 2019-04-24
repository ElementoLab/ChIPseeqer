#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;


my %H1 = ();
my %H2 = ();
my $l = <STDIN>;
print $l;
while (my $l = <STDIN>) {
  
  chomp $l;
  
  my @a = split /\t/, $l, -1;
  
  # remove ID
  my $c = shift @a;
  
  my @d = ();
  my $cnt_pres = 0;
  foreach my $r (@a) {
    if (($r ne "") && ($r ne 'NaN') && ($r ne 'nan')) {
      push @d, $r;
      $cnt_pres ++;
    }
  }
  
  if ($cnt_pres > 0) {
  
    # average 
    my $avg = Sets::average(\@d);
    my $std = Sets::stddev(\@d); 
    
    my @b   = ();
    
    foreach my $r (@a) {
      if (($r ne "") && ($r ne 'NaN') && ($r ne 'nan') && ($r ne "NA")) {
	
	$r = ( $r - $avg );
	if (abs($r) > 0.00000001) {
	  $r = $r / $std;
	} 
	
	push @b, sprintf("%4.3f", $r);
	
      } else {
	
	$r = "";
	push @b, $r;
	
      }
      
    }
    
    print $c . "\t" . join("\t", @b) . "\n";
  } else {
    print "$l\n";
    
  }
}



   
