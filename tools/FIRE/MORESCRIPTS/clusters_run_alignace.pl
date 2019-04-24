#
# 1. get genes that have a given motif
# 2. do a GO analysis on them
#

use Sets;
use Table;
use Getopt::Long;
use Fasta;


use strict;

if (@ARGV == 0) {
  die "USage: perl clusters_run_alignace.pl --clusters=FILE --outfile\n";
}

my $clusters  = undef;
my $fastafile = undef;
my $outfile   = undef;

GetOptions ('clusters=s'        => \$clusters,
	    'fastafile=s'       => \$fastafile,
            'outfile=s'         => \$outfile);



my $fa = Fasta->new;
$fa->setFile($fastafile);

my %SEQ = ();
while (my $a_ref = $fa->nextSeq()) {
    my ($n, $s) = @$a_ref;
    $SEQ{$n} = $s;
}

my $ta = Table->new;
$ta->loadFile($clusters);
my $a_ref = $ta->getArray();
shift @$a_ref;

my %GROUPS        = ();

foreach my $r (@$a_ref) {
  push @{ $GROUPS{ $r->[1] } }, $r->[0];
}

if (defined($outfile)) {
  open OUCH, ">$outfile" or die "Cannot open $outfile.\n";
}

my %H = ();
my $cnt_global = 0;
foreach my $gk (sort {$a <=> $b} (keys(%GROUPS))) {

  my $a_ref_genes = $GROUPS{ $gk };

  my $tmpfile = Sets::getTempFile("TMPFILES/tmp.seq");
  
  open OUT, ">$tmpfile";
  foreach my $g (@$a_ref_genes) {
    print OUT ">$g\n$SEQ{$g}\n\n";
  }	
  close OUT;

  my $out = `alignace2004/AlignACE -i tmp.seq`;
  
  print "Cluster $gk\n\n";
  print "$out\n";
  my $cnt = 0;
  while ($out =~ /MAP\ Score\: ([\d\.]+)/g) {
    if ($1 >= 10) {
      $cnt ++;
    }
  }
  
  $H{ $gk } = $cnt;

  $cnt_global += $cnt;

  
  print "Number of motifs for cluster $gk = $cnt.\n";

  unlink $tmpfile;

}

print "Total number of motifs = $cnt_global.\n";

close OUCH if (defined($outfile));
