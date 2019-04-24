use lib "$ENV{FIREDIR}/SCRIPTS";

use Table;
use Sets;

use strict;

if (@ARGV == 0) {
  die "Usage: perl generate_motif_profiles.pl -fastafile FILE -summaryfile FILE -outfile FILE -rna INT\n";
}

my $progdir     = "$ENV{FIREDIR}/PROGRAMS";
if ((! -e "$progdir/genregexp") && (! -e "$progdir/genregexp.exe")) {
  die "Cannot find genregexp, please check $progdir\n";
}

my $fastafile   = Sets::get_parameter(\@ARGV, "-fastafile");
my $summaryfile = Sets::get_parameter(\@ARGV, "-summaryfile");
my $outfile     = Sets::get_parameter(\@ARGV, "-outfile");

my $rna         = undef;
if (Sets::exist_parameter(\@ARGV, "-rna") == 1) {
  $rna         = Sets::get_parameter(\@ARGV, "-rna");
  print "Setting -rna to $rna.\n";
}

my $noflank     = 0;
if (Sets::exist_parameter(\@ARGV, "-noflank") == 1) {
  $noflank         = Sets::get_parameter(\@ARGV, "-noflank");
  print "Not showing flanking sequence.\n";
}


my $rootdir     = ".";


my $ta = Table->new;
$ta->loadFile($summaryfile);
my $a_ref_mo = $ta->getArray();
my %STAT      = ();
my @MOTIFS    = ();
foreach my $r (@$a_ref_mo) {
  push @MOTIFS, $r->[0];
}


my @HITS = ();
foreach my $re (@MOTIFS) {
  print "Generating profile for $re ... ";
  my $ff = Sets::getTempFile("genregexp");
  my $todo = "$progdir/genregexp -re $re -fastafile $fastafile ";
  
  if ($noflank == 0) {
    $todo .= "-flank 10 ";
  }
  
  if (defined($rna)) {
    $todo .= " -rna $rna ";     
    if ($rna == 1) {
      $todo .= " -reldist 0 ";
    }    
  }
  $todo .= " > $ff"; 
  system($todo);
  $ta->loadFile($ff);
  my $a_ref = $ta->getArray();
  foreach my $r (@$a_ref) {
    my @a = ($re, $r->[0], $r->[1], $r->[2], $r->[3]);
    push @HITS, \@a;
  }
  unlink $ff;
  print "Done.\n";
}

open OUT, ">$outfile" or die "cannot open \"$outfile\".\n";
foreach my $r (@HITS) {
  print OUT join("\t", @$r); print OUT "\n";
}
close OUT;
