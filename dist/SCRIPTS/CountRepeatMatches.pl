#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;


use Getopt::Long;

if (@ARGV == 0) {
  die "Args --table=FILE --samfile=FILE --uniq=INT --name=STR --t2c=INT\n";
}

my $table   = undef;
my $samfile = undef;
my $uniq    = 1;
my $name    = "";
my $t2c     = 0;

GetOptions("table=s"   => \$table,
           "uniq=s"    => \$uniq,
	   "t2c=s"     => \$t2c,
           "name=s"    => \$name,
           "samfile=s" => \$samfile);

open IN, $table or die "Cannot open file $table\n";
my %C = ();
my $numreads = 0;
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  if ($uniq == 0) {
    $C{$a[0]}  = $a[3];
    $numreads += $a[3];
  } else {
    $C{$a[0]}  = 1;
    $numreads ++;
  }

}
close IN;


open IN, $samfile or die "Cannot open $samfile\n";
#$numreads = 0;
my %H = ();
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  
  next if ($a[2] eq "*");

  next if ($a[2] =~ /Simple_repeat/);
  next if ($a[2] =~ /Low_complexity/);

  if (($t2c == 1) && !Sets::hasT2C($a[$#a], $a[9], $a[5]))  {
    next;
  }	


  my @b = split /\-/, $a[2];
  pop @b;
  my $te = join("-", @b);

  $H{$te} += $C{$a[0]};

}
close IN;


print "Repeat\tnumReads$name\tnumReadsLib$name\tfracReads$name\n";
foreach my $k (keys(%H)) {
  my $freq = sprintf("%4.3f", 100 * $H{$k} / $numreads);
  print "$k\t$H{$k}\t$numreads\t$freq\n";
}

