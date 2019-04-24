#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";

use strict;
use Getopt::Long;
use Sets;

my $motif          = undef;
my $expfile        = undef;
my $promcpgislands = undef;
my $promlengths    = "/Users/olivier/PEOPLE/WEIMIN/DESIGN/OLIVIER/promoter_regions_lengths.txt";
my $randomize      = 0;
my $target         = 1;
my $maxnumcpg      = undef;
my $normalize      = 1;

GetOptions("promcpgislands=s" => \$promcpgislands,
           "expfile=s"        => \$expfile,
	   "target=s"         => \$target,
	   "normalize=s"      => \$normalize,
	   "randomize=s"      => \$randomize,
	   "maxnumcpg=s"      => \$maxnumcpg,
	   "motif=s"          => \$motif);

# load target or not
my %lengths = ();
open IN, $promlengths or die "cannot open lens\n";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  $lengths{$a[0]} = $a[1];
}
close IN;


# load target or not
my %targets = ();
open IN, $expfile;
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  $targets{$a[0]} = $a[1];
}
close IN;

# load motif matches
my %MATCHES = ();  # gene x motif arrays of positions
my %MOTIFS  = ();
my $file = Sets::filename($expfile);
open IN, "$expfile\_FIRE/DNA/$file.profiles";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  push @{$MATCHES{$a[1]}{$a[0]}}, $a[2];
  $MOTIFS{$a[0]} = 1; # store motif
}
close IN;




# read CpG island positions in promoters
open IN, $promcpgislands or die "Cannot open $promcpgislands\n";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  my $g = shift @a;
  my $c = shift @a;
  #print "$l\n";
  next if ($c == 0);

  if (defined($maxnumcpg) && ($c > $maxnumcpg)) {
    next;
  }

  #xsprint "Got here target{$g} $targets{$g}\n";
  next if ($targets{$g} == $target);

  next if (!defined($MATCHES{$g}{$motif}));

  if (@{ $MATCHES{$g}{$motif}} > 0) {

    foreach my $p (@{$MATCHES{$g}{$motif}}) {

      if ($randomize == 1) {
	$p = int( 0.5 +  rand($lengths{$g}));	
      }

      foreach my $cpg (@a) {

	my ($i,$j) = split /\-/, $cpg;
	my $l = $j - $i + 1;
	my $relp = undef;

	#print "Comparing $p to [$i,$j]\n";

	if (($i <= $p) && ($p <= $j)) {
	  $relp = Sets::min($p-$i, $j-$p);
	} elsif ($p < $i) {
	  my $d = $i - $p;
	  $relp = - $d;
	} elsif ($p > $j) {
	  my $d = $p - $j;
	  $relp = - $d;
	}
	
	if ($normalize == 1) {
	  $relp /= $l;
	}

	print "$relp\n";

      }
    }
  }

}
close IN;


