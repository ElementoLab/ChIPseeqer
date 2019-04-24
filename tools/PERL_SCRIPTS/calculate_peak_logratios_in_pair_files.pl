#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Statistics::Distributions;
use Sets;
#use strict;

if (@ARGV == 0) {
  die "Args: suffix (will add _532.pair)\n";
}

my $suf = $ARGV[0];

srand(1234);

#
# load 532 data
#
my $f532 = "$suf\_532.pair";
if (! -e $f532) {
  $f532 .= ".txt";
}
open IN, $f532 or die "Cannot open $f532\n";
my $l = <IN>;
my $l = <IN>;

my %PROMS532 = ();

while (my $l = <IN>) {
  next if ($l =~ /^\#/);
  chomp $l;
  my @a = split /\t/, $l, -1;
  next if ($a[2] =~ /RANDOM/);
  if ($a[9] eq "") {
    die "At $a[3], no intensity.\n";
  }
  my @probes = ($a[3], $a[9]);
  push @{ $PROMS532{$a[2]} }, \@probes;
}
close IN;

#
# load 635 data
#
my $f635 = "$suf\_635.pair";
if (! -e $f635) {
  $f635 .= ".txt";
}
open IN, $f635 or die "Cannot open $f635\n";;
my $l = <IN>;
my $l = <IN>;

my %PROMS635 = ();

while (my $l = <IN>) {
  next if ($l =~ /^\#/);
  chomp $l;
  my @a = split /\t/, $l, -1;
  next if ($a[2] =~ /RANDOM/);
  my @probes = ($a[3], $a[9]);
  if ($a[9] eq "") {
    die "At $a[3], no intensity.\n";
  }	
  push @{ $PROMS635{$a[2]} }, \@probes;
}
close IN;


#
# calculate ratio profile
#

my @intensities = ();
my %cprobes     = ();  # combine probes
my @nprobes     = ();

foreach my $p (keys(%PROMS532)) {

  my @a_532 = sort sortProbes @{$PROMS532{$p}};
  my @a_635 = sort sortProbes @{$PROMS635{$p}};
  
  my $n = @a_532;

  my @v = ();
  for (my $i=0; $i<$n; $i++) {
    if ($a_532[$i]->[0] ne $a_635[$i]->[0]) {
      die "Problem with data for $p\n";
    }

    my $l = $a_635[$i]->[1] / $a_532[$i]->[1];
    push @v, $l;
    push @intensities, $l;
  }
  
  $cprobes{$p} = \@v;
  push @nprobes, scalar(@v);
}

my $ng = @nprobes;
my $ni = @intensities;



open OUT, ">$suf\_peaks.tmp" or die "Cannot open $suf\_peaks.tmp";
foreach my $p (keys(%cprobes)) {
  #next if ($p ne "HSAP0406S00022271");

  my $np = scalar(@{$cprobes{$p}});

  next if ($np < 10);
  

  my $a_ref_s = &getSmoothProbeProfile($cprobes{$p});
  my $m       = Sets::maxInArray($a_ref_s);


  #my $z   = ($m - $avg) / $std;
  #my $uprob = Statistics::Distributions::uprob ($z);


  my @a_tmp = ($p, sprintf("%4.3f", $m), $np); # sprintf("%4.3f", $z), $uprob);
  
  push @pvall, \@a_tmp;
 
  print OUT join("\t", @a_tmp) . "\n";
  
 #push @pv,    $uprob;

}
close OUT;
print "$suf\_peaks.tmp created";



sub sortProbes {
  my ($p1) = $a->[0] =~ /^CHR\d+PR{0,1}0*(\d+?)$/;
  my ($p2) = $b->[0] =~ /^CHR\d+PR{0,1}0*(\d+?)$/;
  #print "$p1 <=> $p2\n";
  return $p1 <=> $p2;
}

sub getSmoothProbeProfile {
  my ($a_ref_p) = @_;

  my $n = @$a_ref_p;
  my @a = ();
  for (my $i=0; $i<$n-2; $i++) {
    my $sum = 0;
    for (my $j=$i; $j<$i+3; $j++) {
      $sum += $a_ref_p->[$j];
    } 
    $sum /= 3;
    $a[$i] = $sum;
  }
  return \@a;
}


sub getGEVpvalues {
  my ($nullfile, $maxfile) = @_;
  
  my $txt = "
library(evd)
library(ismev)
nuld <- read.table(\"$nullfile\")
v    <- gev.fit(nuld[,1], show=F)\$mle
maxd <- read.table(\"$maxfile\")
outd <- pgev(maxd[,2], v[1], v[2], v[3], lower.tail=F)
outm <- p.adjust(outd, method=\"BH\")
write.table(cbind(outd,outm), file=\"$maxfile.out\", sep=\"\\t\", row.names=T, col.names=NA, quote=F)
";
  #print "$txt\n";
  my $ft = Sets::getTempFile("/tmp/Rscript");
  Sets::writeText($txt, $ft);
  if (! -e $ft) {
    die "Problem creating $ft\n";
  }
  my $todo = "/usr/bin/R CMD BATCH $ft";
  #print STDERR "$todo\n";
  system($todo) == 0 or die "Cannot exec R script ?\n";
  
  my @b = ();
  open INCH, "$maxfile.out" or die "Cannot open $maxfile.out ... verify R script\n";
  my $l = <INCH>;
  while (my $l = <INCH>) {
    chomp $l;
    my @a = split /\t/, $l, -1;
    shift @a;
    push @b, \@a;    
  }
  close INCH;

  return \@b;

}

