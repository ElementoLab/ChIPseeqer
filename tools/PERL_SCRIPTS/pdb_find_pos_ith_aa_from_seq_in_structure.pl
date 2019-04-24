#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;
use Fasta;
use ClustalW;

use Getopt::Long;

if (@ARGV == 0) {
  die "Args --NM=STR --pos=INT --pdbfa=FILE [ --protfile=FILE ] \n";
}

my $protfile = "$ENV{SNVSEEQERDIR}/REFDATA/human.protein.clean.faa";
my $NM       = undef;
my $pos      = undef;
my $pdbfa    = undef;

GetOptions("NM=s"       => \$NM,
	   "protfile=s" => \$protfile,
           "pdbfa=s"    => \$pdbfa,
	   "pos=s"      => \$pos);


my $np = undef;
open IN, "$ENV{SNVSEEQERDIR}/REFDATA/refLinkrefGene.txt.26Dec2009" or die "Cannot open file \n";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  if ($a[0] eq $NM) {
    $np = $a[1];
  }
}
close IN;

if (!defined($np)) {
  die "Cannot find $NM\n";
}

# get protein seq (NP)
my @a = `fastacmd -d $protfile -s $np`;
shift @a;
chomp @a;
my $txt = join("", @a);
#print "$txt\n";

# 
my $fa = Fasta->new;
$fa->setFile($pdbfa);
my $a_ref = $fa->nextSeq();
my ($n, $seq) = @$a_ref;
$seq =~ s/\-//g;
$seq =~ s/\*//g;

# store
my $tmpfile = Sets::getTempFile("/tmp/NPprot");
open OUT, ">$tmpfile";
print OUT ">$np\n$txt\n";
print OUT ">MOD\n$seq\n";
close OUT;

# align
my $todo = "clustalw $tmpfile > /dev/null";
system($todo) == 0 or die "Cannot exec $todo\n";

# 
my $cl = ClustalW->new;

$cl->setFile("$tmpfile.aln");

my $a_ref_aln = $cl->getSeqsWithNames();

my $s1 = $a_ref_aln->[0]->[1];
my $s2 = $a_ref_aln->[1]->[1];

#print "$s1\n";
#print "$s2\n";

my @a1 = split //, $s1;
my @a2 = split //, $s2;

while ($a1[0] eq '-') {
  shift @a1;
  shift @a2;
}

#print join("", @a1) . "\n";

my $i1 = 0;
my $i2 = 0;
for (my $j=0; $j<@a1; $j++) {
  
  #print "$a1[$j]$a2[$j]\n";

  if ($a1[$j] ne '-') {
    $i1++;
  }
  if ($a2[$j] ne '-') {
    $i2++;
  }

  if ($i1 == $pos) {
    print "$i2 ($a2[$j])\n";
    last;
  }


}





#my @b = split //, $txt;
#print $b[$ARGV[2]-1] . "\n";

unlink $tmpfile;
