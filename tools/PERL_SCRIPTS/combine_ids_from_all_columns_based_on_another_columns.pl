BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my %H = ();
foreach my $r (@$a_ref) {  
  #next if ($r->[$ARGV[1]] eq "");
  #next if ($r->[$ARGV[2]] eq "");

  my $n = shift @$r;

  for (my $i=0; $i<@$r; $i++) {
    #print "push $r->[$i]\n";
    push @{ $H{ $n }->[$i] }, $r->[$i] if (($r->[$i] ne "") && (!Sets::in_array($r->[$i], @{ $H{ $n }->[$i] })));

  }
}


foreach my $k (keys(%H)) {
  my $r = $H{$k};
   
  print "$k";
  for (my $i=0; $i<@$r; $i++) {
    print "\t" . join("/", @{$H{$k}->[$i]});
  }	
  print "\n";
}



