BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use FileHandle;
use Table;
use Sets;

my $ta = Table->new;

my $a_ref = Sets::readSet($ARGV[0]);
shift @ARGV;

my @FH = ();
foreach my $f (@ARGV) {
  
  $ta->loadFile($f);
  my $h_ref = $ta->getIndexKV(0,1);
  
  foreach my $s (@$a_ref) {
    push @{ $H{ $s } }, $h_ref->{$s};
  }

}

foreach my $s (sort { $a <=> $b } @$a_ref) {
  print "$s\t" . join("\t",  @{ $H{ $s } }) . "\n"; 
}


