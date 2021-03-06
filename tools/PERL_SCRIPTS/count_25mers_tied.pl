BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Fasta;
use Sets;
use DB_File;
use Fcntl;    # For the constants O_RDWR and O_CREAT
use strict;

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);

my %C = ();

tie (%C, 'DB_File', "counts.dat", O_RDWR|O_CREAT, 0666, $DB_BTREE);

while (my $a_ref = $fa->nextSeq()) {
    my ($n, $s) = @$a_ref;
    my $l = length($s);
    for (my $i=0; $i<$l-25; $i++) {
      
      my $ss = substr($s, $i, 25);
      $C{ $ss } ++;
      my $sc = Sets::getComplement($ss);
      $C{ $sc } ++;
      
    }
}


my $a_ref = Sets::hash_order(\%C);

my %CC = ();
foreach my $r (@$a_ref) {
  $CC{ $C{$r} } ++;
}

foreach my $r (sort { $a <=> $b } (keys(%CC))) {
  print "$r\t$CC{$r}\n";
}

untie %C;

