#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;
use strict;

if (@ARGV == 0) {
  die "Args: fastafile\n";
}

my $a_ref_j = Sets::getFiles("$ENV{HOME}/PROGRAMS/MYSCANACE/JASPAR/*.jaspar");

my $n = @$a_ref_j;

print STDERR "Found $n matrices\n";

my %MATRIX = ();
my @MOT    = ();
foreach my $m (@$a_ref_j) {

  print STDERR "$m\n";
  my $ff   = Sets::filename($m);
  $ff =~ s/\.jaspar//;
  push @MOT, $ff;

  $m =~ s/\(/\\\(/g;
  $m =~ s/\)/\\\)/g;

  my $todo = "$ENV{HOME}/PROGRAMS/MYSCANACE/MyScanACE -z $ARGV[0] -g 0.52 -c 1 -g 0.5 -allowmatchoverlap 0 -p 1 -output density -header 0 ";
  $todo .= " -j $m ";
  #print "$todo\n";
  my $txt  = `$todo`;
  
  my @a    = split /\n/, $txt;
  foreach my $l (@a) {
    my @b = split /\t/, $l;
    $MATRIX{$b[0]}{$ff} = $b[1];  #  gene x motif = density
  }	
}


my @genes = keys(%MATRIX);

# determine which of the motifs have [ 5% 95% ] match 
my %EXC = ();
foreach my $s (@MOT) {
  my $numgeneswithmatch = 0;
  foreach my $g (@genes) {
    if ($MATRIX{$g}{$s} > 0) {
      $numgeneswithmatch ++;
    }
  }

  my $fa = $numgeneswithmatch / scalar(@genes);

  if (($numgeneswithmatch < 100) || ($fa > 0.9)) {
    $EXC{$s} = 1;
  }

}

print "G";
foreach my $m (@MOT) {
  next if (defined($EXC{$m}));
  print "\t$m";
}
print "\n";

foreach my $g (@genes) {
  print "$g";
  foreach my $s (@MOT) {
    next if (defined($EXC{$s}));
    print "\t" . $MATRIX{$g}{$s};
  }
  print "\n";
}
