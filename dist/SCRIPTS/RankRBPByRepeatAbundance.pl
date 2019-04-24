#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;

use Getopt::Long;

if (@ARGV == 0) {
  die "Args --f1=FILE --f2=FILE\n";
}
my $files   = undef;
my $repeats = undef;
my $RBPs    = undef;
my $exclRBP = undef;
my $families = 0;

GetOptions("files=s"   => \$files,
	   "RBPs=s"    => \$RBPs,
	   "exclRBP=s" => \$exclRBP,
	   "families=s"=> \$families,
           "repeats=s" => \$repeats);


# make list of repeats
my @a_rep   = ();

if ($families == 0) {

  my $replist = "/panda_scratch_miro/ole2001/STAU2/repeatnames.txt";  
  open IN, $replist or die "Cannot open file $replist";
  while (my $l = <IN>) {
    chomp $l;
    my @a = split /\t/, $l, -1;
    if ($l =~ /$repeats/) {
      push @a_rep, $l;
    }
  }
  close IN;

} else {
  
  open IN, "$ENV{CHIPSEEQERDIR}/DATA/RepMask3.2.7_annotation_families.txt" or die "Cannot open file";
  while (my $l = <IN>) {
    chomp $l;
    my @a = split /\t/, $l, -1;
    next if ($a[0] =~ /\?$/);
    push @a_rep, $a[0];
  }
  close IN;
  
}

print STDERR "Found " . scalar(@a_rep) . " matching repeat families\n";

my %RANKS = ();
my %RBP   = (); 
foreach my $rep (@a_rep) {
  
  my $todo = "perl /home/ole2001/PROGRAMS/ChIPseeqer-1.0/SCRIPTS/SummarizeRepeatFamilyAbundanceInLibraries.pl --files=\"$files\" --type=\"$rep\" ";
  if ($families == 1) {
    $todo .= " --fam=1 ";
  }
  $todo .= " | sort_column_inv.pl 1 | columns.pl 1 0";
  
  my $txt = `$todo`;
  my @a_lines = split /\n/, $txt;
  for (my $i=0; $i<@a_lines; $i++) {
    my @a = split /\t/, $a_lines[$i];
    if ($a[1] =~ /$RBPs/) {

      next if (defined($exclRBP) && ($a[1] eq $exclRBP));

      $RANKS{$rep}{$a[1]} = $i+1;
      $RBP  {$a[1]} = 1;
    }
  }	
  
}

foreach my $r (keys(%RBP)) {
  print "\t" . $r;
}
print "\n";

foreach my $rep (@a_rep) {
  print "$rep";
  foreach my $r (keys(%RBP)) {
    print "\t" . $RANKS{$rep}{$r};
  }
  print "\n";
}



