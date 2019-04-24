use Getopt::Long;

my $utr3file  = "$ENV{SNVSEEQERDIR}/REFDATA/refGene.txt.25Nov2009.utr3seq";
my $seedfile  = "$ENV{SNVSEEQERDIR}/REFDATA/mirnas-29DEC2009.seeds";
my $stopafter = undef;
 
GetOptions("utr3file=s"  => \$utr3file,
	   "stopafter=s" => \$stopafter,	
           "seedfile=s"  => \$seedfile);


my %MATRIX = ();

open IN, $seedfile;
my @SEEDS = ();
my $num = 0;
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;

  my $todo = "$ENV{SNVSEEQERDIR}/genregexp -re $a[0] -fastafile $utr3file -rna 1 -format density";
  my $seedname = "$a[0]-$a[2]";
  my $txt  = `$todo`;
  my @lines = split /\n/, $txt;
  foreach my $l (@lines) {
    my @b = split /\t/, $l;
    $MATRIX{$b[0]}{$seedname} = $b[1];
  }
 
  push @SEEDS, $seedname;

  $num ++;
  if (defined($stopafter) && ($num == $stopafter)) {
    last;
  }
}
close IN;

print "GENE";
foreach my $s (@SEEDS) {
  print "\t$s";
}
print "\n";

foreach my $gene (keys(%MATRIX))  {

  print "$gene";
  foreach my $s (@SEEDS) {
    print "\t$MATRIX{$gene}{$s}";
  }
  print "\n";
}
  
