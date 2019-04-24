BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;
use strict;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my $h_ref = $ta->getIndexShifted();

my $p = $h_ref->{$ARGV[1]};

shift @$a_ref;

my %H = ();
foreach my $r (@$a_ref) {
  #print ">$r->[0]\n";
  #print join("\t", @$p); print "\n";
  #print join("\t", @{ $h_ref->{$r->[0]} }); print "\n";
  
  #if (Sets::pearson($p, $h_ref->{$r->[0]}) > 0.8) {
   # print "$r->[0]\n";
    #my @p1 = @{$r}[1..10];
    #my @p2 = @{$r}[11..18];
     
    my @p1 = @{$r}[1..3]; 
    my @p2 = @{$r}[4..6];

    #print join(" ", @p1); print "\n";
    #print join(" ", @p2); print "\n";
    
    $H{ $r->[0] } = Sets::pearson(\@p1, \@p2);
  #}
}

my $o_ref = Sets::hash_order(\%H);
foreach my $k (reverse(@$o_ref)) {
  print "$k\t" . sprintf("%4.3f", $H{$k}) . "\n";
}

