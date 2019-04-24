BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Fasta;
use Table;

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);
my $a_ref = $fa->nextSeq();
my ($n, $s) = @$a_ref;



my $ta = Table->new;
$ta->loadFile($ARGV[1]);
my $a_ref = $ta->getArray();

foreach my $r (@$a_ref) {
  my $l = $r->[3] - $r->[2] + 1;
  my $news = 'N' x $l;
  substr($s, $r->[2], $l) = $news;
}

my @a = split /N+/, $s;

my $cnt = 1;
foreach my $s (@a) {
  if (length($s) > 10) {
    print ">SEQ$cnt\n$s\n\n";
    $cnt ++;
  }
}


