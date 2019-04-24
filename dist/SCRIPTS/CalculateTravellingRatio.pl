use lib "$ENV{CHIPSEEQERDIR}/SCRIPTS";
use strict;
use Sets;
use Getopt::Long;

my $lenu     = 30;
my $lend     = 300;
my $peakfile = undef;
my $chipdir  = undef;
my $verbose  = 0;
my $fraglen  = 0;

GetOptions("lenu=s"     => \$lenu,
	   "lend=s"     => \$lend,
 	   "fraglen=s"  => \$fraglen,
	   "verbose=s"  => \$verbose,
	   "chipdir=s"  => \$chipdir,
	   "peakfile=s" => \$peakfile);


my $todo0 = "ChIPseeqerSummary --targets=$peakfile --lenu=$lenu --lend=$lend";
print "$todo0\n" if ($verbose == 1);
system($todo0) == 0 or die "cannot exec $todo0\n";
my $outNM = "$peakfile.RefGene.NM.txt";
open IN, $outNM or die "Cannot open $outNM\n";
my %BOUNDTS = ();
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  $BOUNDTS{$a[0]} = 1;
}
close IN;


#
# get prox intervals
#
my $tmpfile1 = Sets::getTmpFile("prox");
my $todo1 = "perl $ENV{CHIPSEEQERDIR}/SCRIPTS/extract_upstream_sequence_coordinates_from_annotation.pl --annotation=$ENV{CHIPSEEQERDIR}/DATA/refGene.txt.25Nov2009.oneperTSS --checkmaxlen=0 --lengthU=$lenu --lengthD=$lend > $tmpfile1";
print "$todo1\n" if ($verbose == 1);
system($todo1) == 0 or die "cannot exec $todo1\n";
open IN, $tmpfile1;
open OUT, ">$tmpfile1.bound";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l;
  if ($BOUNDTS{$a[0]}) {
    print OUT "$a[0]\t$a[1]\t$a[2]\t$a[3]\n";
  }
}
close IN;
close OUT;
if (-e "$tmpfile1.bound") {
  print "$tmpfile1.bound created.\n";
}


# get avg read count
$todo1 = "$ENV{CHIPSEEQERDIR}/ChIPgetIntervalReadCounts -hasid 1 -intervals $tmpfile1.bound -fraglen $fraglen -normalize 0 -chipdir $chipdir -output avg -chrdata $ENV{CHIPSEEQERDIR}/DATA/hg18.chrdata > $peakfile.avgprox";
print "$todo1\n" if ($verbose == 1);
system($todo1) == 0 or die "Cann exec $todo1\n";



#
# get body intervals
#
my $tmpfile2 = Sets::getTmpFile("body");
my $todo2 = "perl  $ENV{CHIPSEEQERDIR}/SCRIPTS/extract_upstream_sequence_coordinates_from_annotation.pl --annotation=$ENV{CHIPSEEQERDIR}/DATA/refGene.txt.25Nov2009.oneperTSS --checkmaxlen=0 --lengthU=-$lend --lengthD=TES > $tmpfile2";
print "$todo2\n" if ($verbose == 1);
system($todo2) == 0 or die "cannot exec $todo2\n";
open IN, $tmpfile2;
open OUT, ">$tmpfile2.bound";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l;
  if ($BOUNDTS{$a[0]}) {
    print OUT "$a[0]\t$a[1]\t$a[2]\t$a[3]\n";
  }
}
close IN;
close OUT;
if (-e $tmpfile2) {
  print "$tmpfile2 created.\n";
}

# get avg read count
$todo2 = "$ENV{CHIPSEEQERDIR}/ChIPgetIntervalReadCounts -hasid 1 -intervals $tmpfile2.bound -fraglen $fraglen -normalize 0 -chipdir $chipdir -output avg -chrdata $ENV{CHIPSEEQERDIR}/DATA/hg18.chrdata > $peakfile.avgbody";
print "$todo2\n" if ($verbose == 1);
system($todo2) == 0 or die "Cann exec $todo2\n";

# ratio

my %rbody = ();
open IN, "$peakfile.avgbody";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  $rbody{$a[0]} = $a[4];
}
close IN;

open IN, "$peakfile.avgprox";
open OUT, ">$peakfile.TR";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  if (defined($rbody{$a[0]})) {
    my $ra = Sets::max($a[4],0.001) / Sets::max($rbody{$a[0]}, 0.001);
    print OUT "$a[0]\t$a[4]\t$rbody{$a[0]}\t$ra\n";
  }
}
close IN;
close OUT;
print "$peakfile.TR created\n";
