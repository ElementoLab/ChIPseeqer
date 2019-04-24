use lib qw(/home/elemento/PERL_MODULES);
use Sets;
use strict;

open IN, $ARGV[0];
my $g = undef;
my %H = ();

while (my $l = <IN>) {
  chomp $l;
  

  if ($l =~ /\-\-\>\t(.+?)$/) {
    $g = $1;
  } else {
    
    my @a = split /\t/, $l, -1;
    
    if ($a[1] >= $ARGV[1]) {
      
      my @b = split /\//, $a[2];
      
      foreach my $r (@b) {
	push @{ $H{ $r } }, $g if (!Sets::in_array($g, @{ $H{ $r } }));
      }
      
    }
    
  }  
}

foreach my $r (keys(%H)) {
  print "$r\t"; print join("\t", @{ $H{ $r } }); print "\n";
}
