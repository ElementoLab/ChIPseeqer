BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Fasta;
use Getopt::Long;
use Table;
use strict;

my $fasta1    = undef;
my $fasta2    = undef;
my $outdir    = undef;
my $orthologs = undef;

GetOptions (
	    'fasta1=s'    => \$fasta1,
	    'fasta2=s'    => \$fasta2,
	    'outdir=s'    => \$outdir,
	    'orthologs=s' => \$orthologs
	   );

my $fa = Fasta->new;
$fa->setFile($fasta1);

my %SEQ1 = ();
while (my $a_ref = $fa->nextSeq()) {
  my ($n, $s) = @$a_ref;
  $SEQ1{$n} = $s;
}

$fa->dispose();

$fa = Fasta->new;
$fa->setFile($fasta2);

my %SEQ2 = ();
while (my $a_ref = $fa->nextSeq()) {
  my ($n, $s) = @$a_ref;
  $SEQ2{$n} = $s;    
}



my $ta = Table->new;
$ta->loadFile($orthologs);
my $a_ref = $ta->getArray();

my $f1 = undef;
my $f2 = undef;

if (defined($outdir)) {

  my $fn1 = Sets::filename($fasta1);
  $f1 = "$outdir/$fn1.fc1";
  my $fn2 = Sets::filename($fasta2);
  $f2 = "$outdir/$fn2.fc2";

} else {

  $f1 = "$fasta1.fc1";
  $f2 = "$fasta2.fc2";

}


open OUT1, ">$f1";
open OUT2, ">$f2";

foreach my $r (@$a_ref) {

  my @a = @$r;
  if (!defined($a[1])) {
    $a[1] = $a[0];
  }

  

  print OUT1 ">$a[0]\n$SEQ1{$a[0]}\n\n";

  

  print OUT2 ">$a[1]\n$SEQ2{$a[1]}\n\n";

}

close OUT1;
close OUT2;

