BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;
use Table;
use File::Copy;
use strict;

my $ta = Table->new;
#$ta->setLimit(1);
$ta->loadFile($ARGV[0]) or "Cannot load matrix1 file\n";
my $a_ref1 = $ta->getArray();


my $a_ref2 = Sets::readSet($ARGV[1]);  # sorted
my @cols   = @{ $a_ref1->[0] }; shift @cols;


foreach my $r (@$a_ref1) {

  
  my $i      = 0;
  my %H      = ();
  foreach my $c (@$a_ref2) {
    push @{ $H{ $c } }, $i; 
    $i ++;
  }

  #foreach my $k (keys(%H)) {
  #  print "$k => " . join("=\t=", @{ $H{$k} }) . "\n";
  #}

  my $n = shift @$r;
  print "$n\t";

  my @newcol = ();
  for (my $i=0; $i<@$r; $i++) {
    my $reali = shift @{ $H{ $cols[$i] } };
    $newcol[$reali] = $r->[$i];
    #print "newcol[$reali] = $r->[$i];\n";
  }

  
  print join("\t", @newcol) . "\n";
  
}
  


