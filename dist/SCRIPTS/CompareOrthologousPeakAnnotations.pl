use Getopt::Long;
use strict;

my $sp1annot  = undef;
my $sp2annot  = undef;
my $orthology = undef;
my $revorth   = 0;
my $allts     = undef;

GetOptions("sp1annot=s"    => \$sp1annot,
           "sp2annot=s"    => \$sp2annot,
	   "orthology=s"   => \$orthology,
	   "allts=s"       => \$allts,
	   "revorth=s"     => \$revorth);

# fraction of conserved tgt mouse genes w peaks in sp2 / conserved tgt mouse genes
# human tgt genes with mouse cons vs human genes with mouse cons

# num ts
my $numtssp1 = 0;
open IN, $allts;
while (my $l = <IN>) {
  $numtssp1 ++;
}
close IN;



my %ORTH = ();
my %ORTHREV = ();

my $numconstssp2 = 0;

open IN, $orthology;
my $numtssp1withorth = 0;
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;

  next if ($a[1] =~ /reciprocal/);
  next if ($a[1] =~ /found/);

  if ($revorth == 1) {
    $ORTH{$a[1]} = $a[0];
    $ORTHREV{$a[0]} = $a[1];
    $numtssp1withorth++;
  } else {
    $ORTH{$a[0]} = $a[1];
  }
  
  $numconstssp2++;
  
}
close IN;

my $numconstssp2_tgt = 0;
# open peak annot for sp2
open IN, $sp2annot;
my $h2 = <IN>; chomp $h2;

my %ANNOT2 = ();
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  $ANNOT2{$a[0]} = \@a;
  if (defined($ORTHREV{$a[0]})) {
    $numconstssp2_tgt++;
  }
}
close IN;



# open peak annot for sp1
open IN, $sp1annot;
my $h1 = <IN>; chomp $h1;
print "$h1\t$h2\n";

my $numtgtwithorth = 0;
my $numtgtwithorthcons = 0;

while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  
  my $og = $ORTH{$a[0]};
  if (defined($og)) {
    $numtgtwithorth ++;
  }

  print "$l";
  if (defined($ANNOT2{$og})) {
    print "\t" . join("\t", @{$ANNOT2{ $og } });
    $numtgtwithorthcons ++;
  } elsif (!defined($og)) {
    print "\tNo clear ortholog";
  } else {
    print "\tOrtholog but not peak";
  }
  print "\n";

}
close IN;

my $ra1 = sprintf("%4.3f", $numtgtwithorthcons/$numtgtwithorth);
my $ra2 = sprintf("%4.3f",  $numconstssp2_tgt/$numconstssp2);
print STDERR "$numtgtwithorthcons/$numtgtwithorth = $ra1 .. vs .. $numconstssp2_tgt/$numconstssp2 = $ra2\n";



