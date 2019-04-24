#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;

use Getopt::Long;

if (@ARGV == 0) {
  die "Args: --name=FILE --desc=FILE --gtf=FILE --visibility=INT\n";
}
my $name       = "OETrack";
my $desc       = "OETrack";
my $visibility = 3;
my $gtf        = undef;

GetOptions("gtf=s"        => \$gtf,
	   "name=s"       => \$name,
           "desc=s"       => \$desc,
	   "visibility=s" => \$visibility);




open OUT, ">$gtf.tmp";
open IN, $gtf or die "Cannot open $gtf\n";
print OUT "track name=$name description=\"$desc\" visibility=$visibility\n";
while (my $l = <IN>) {
  #print $l;
  print OUT $l;
}
close IN;
close OUT;
system("cp $gtf $gtf.bak") == 0 or die "Cnnot vreate bckip\n";
system("mv $gtf.tmp $gtf") == 0 or die "Canot creae backup\n";

