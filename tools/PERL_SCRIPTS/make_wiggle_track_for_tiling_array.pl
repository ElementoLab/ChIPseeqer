#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Statistics::Distributions;
use Sets;
use strict;

if (@ARGV == 0) {
  die "Args: suffix translate (will add _532.pair)\n";
}

my $suf = $ARGV[0];
my $wsm = 5;

srand(1234);


#
# read oligo translation file
#
my %TRAN = ();
open IN, $ARGV[1];
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  $TRAN{$a[0]} = $a[1];
}
close IN;



#
# load 532 data
#
my $f532 = "$suf\_532.pair";
if (! -e $f532) {
  $f532 .= ".txt";
}
open IN, $f532 or die "Cannot open $f532\n";
my $l = <IN>;

my %PROMS532 = ();

while (my $l = <IN>) {
  next if ($l =~ /^\#/);
  chomp $l;
  my @a = split /\t/, $l, -1;
  next if ($a[2] =~ /RANDOM/);
  next if ($a[0] eq "IMAGE_ID");
  next if ($a[1] eq "HELP");

  if (defined($TRAN{$a[3]})) {
    $a[3] = $TRAN{$a[3]};
  } else {
    next;
  }

  if ($a[9] eq "") {
    die "At $a[3], no intensity.\n";
  }
  my @probes = ($a[3], $a[9]);
  push @{ $PROMS532{$a[2]} }, \@probes;
}
close IN;

print STDERR "Done reading $f532.\n";

#
# load 635 data
#
my $f635 = "$suf\_635.pair";
if (! -e $f635) {
  $f635 .= ".txt";
}
open IN, $f635 or die "Cannot open $f635\n";;
my $l = <IN>;


my %PROMS635 = ();

while (my $l = <IN>) {
  next if ($l =~ /^\#/);
  chomp $l;
  my @a = split /\t/, $l, -1;
  next if ($a[2] =~ /RANDOM/);
  next if ($a[0] eq "IMAGE_ID");
  next if ($a[1] eq "HELP");
  
  if (defined($TRAN{$a[3]})) {
    $a[3] = $TRAN{$a[3]};
  } else {
    next;
  }


  my @probes = ($a[3], $a[9]);
  if ($a[9] eq "") {
    die "At $a[3], no intensity.\n";
  }	
  push @{ $PROMS635{$a[2]} }, \@probes;
}
close IN;

print STDERR "Done reading $f635.\n";

#
# calculate ratio profile
#

my @intensities = ();
my %cprobes     = ();  # combined probes
my %sprobes     = ();  # same, shuffled
my @probenames  = ();
my @loci        = ();
my @nprobes     = ();


print "track type=wiggle_0 name=\"EZH2_raw\" description=\"EZH2_raw_logratios\" visibility=full autoScale=off\n";


my $prevp = "XXXXXXXXX";
foreach my $p (keys(%PROMS532)) {

  #print "Locus $p\n";

  my @a_532 = sort sortProbes @{$PROMS532{$p}};
  my @a_635 = sort sortProbes @{$PROMS635{$p}};
  
  my $n = @a_532;

  my @v = ();
  for (my $i=0; $i<$n; $i++) {
  
    if ($a_532[$i]->[0] ne $a_635[$i]->[0]) {
      die "Problem with data for $p\n";
    }

    next if ($a_532[$i]->[0] eq $prevp);

    #print " $a_635[$i]->[1] / $a_532[$i]->[1]\n";

    my $l = log( $a_635[$i]->[1] / $a_532[$i]->[1] );

    # push vector of name,logratio

    my ($c, $p) = $a_532[$i]->[0] =~ /^(CHR\d+)PR{0,1}0*(\d+?)$/;

    $c = lc($c);
    $p = $p + 25;
    my $p2 = $p + 1;

    print "$c\t$p\t$p2\t$l\n";

    $prevp =  $a_532[$i]->[0];

  }

}




sub sortProbes {
  my ($p1) = $a->[0] =~ /^CHR\d+PR{0,1}0*(\d+?)$/;
  my ($p2) = $b->[0] =~ /^CHR\d+PR{0,1}0*(\d+?)$/;
  #print "$p1 <=> $p2\n";
  return $p1 <=> $p2;
}

sub getSmoothProbeProfile {
  my ($a_ref_p, $w) = @_;

  my $n = @$a_ref_p;
  my @a = ();
  for (my $i=0; $i<$n-$w+1; $i++) {
    my $sum = 0;
    for (my $j=$i; $j<$i+$w; $j++) {
      if (($a_ref_p->[$j]->[1] eq "") || !defined($a_ref_p->[$j]->[1])) {
	die "Problem with smoothing.\n";
      }
      $sum += $a_ref_p->[$j]->[1];
      
    } 
    $sum /= $w;
    
    # name of middle probe, smoothed avg
    $a[$i] = [ $a_ref_p->[$i+ int(0.5+$w/2) ]->[0], $sum ];
  }
  return \@a;
}


sub getGEVpvalues {
  my ($datafile) = @_;
  
  my $txt = "
library(evd)
library(ismev)
da <- read.table(\"$datafile\")
v    <- gev.fit(da[,4], show=F)\$mle
outd <- pgev(da[,3], v[1], v[2], v[3], lower.tail=F)
outm <- p.adjust(outd, method=\"BH\")
write.table(cbind(da,outd,outm), file=\"$datafile.out\", sep=\"\\t\", row.names=T, col.names=NA, quote=F)
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
  open INCH, "$datafile.out" or die "Cannot open $datafile.out ... verify R script\n";
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

