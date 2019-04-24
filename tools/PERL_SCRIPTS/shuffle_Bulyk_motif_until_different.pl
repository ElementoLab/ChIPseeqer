#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;
use Getopt::Long;

my $cormax = 1/3;
my $motif  = undef;

if (@ARGV == 0) {
  die "Args --motif=FILE --cormax=FILE\n";
}
GetOptions("motif=s"  => \$motif,
           "cormax=s" => \$cormax);


system("cat $motif");
my $tmpfile = Sets::getTempFile("/tmp/motmot");

my $cor = 1;
my $cnt = 0;
do {
  my $todo    = "perl $ENV{HOME}/PERL_SCRIPTS/shuffle_Bulyk_motif.pl $motif > $tmpfile";
  system($todo) == 0 or die "Cannot exec $todo\n";

  #system("cat $tmpfile");
  $todo = "$ENV{HOME}/PROGRAMS/COMPAREACE/MyCompareAce $motif $tmpfile -bulyk ";
  my $txt = `$todo`; 
  #print "$txt\n";
  $txt =~ s/\n//g;
  $cor = $txt;
  print STDERR "Generated randomized matrix with cor = $cor\n";
  $cnt ++;
} while (($cor > $cormax) && ($cnt < 10));

print STDERR "Found a shuffled motif with CompareACE score < $cormax = $cor\n"; 

system("cp $tmpfile $motif.shuffled");
print STDERR "Created $motif.shuffled\n";
system("cat $tmpfile");
unlink $tmpfile;
