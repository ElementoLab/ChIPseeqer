#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;

use Table;


my $ta = Table->new;


# open combines matrix
open IN, $ARGV[0];

my $l = <IN>;
my @motifs = split /\t/, $l, -1;
shift @motifs;

# fix names
foreach my $m (@motifs) {
  $m =~ s/.matches.*$//;
}

# cycle thru genesets
while (my $l = <IN>) {

  chomp $l;
  my @a = split /\t/, $l, -1;

  my $gs = shift @a;
  #next if ($gs !~ /Pf/);
  $gs =~ s/\ .*$//;
  my @mo = ();
  for (my $i=0; $i<@a; $i++) {
    if ($a[$i] > 0.301) {
      push @mo, $motifs[$i];
    }
  }

  print STDERR "Found " . join(" ", @mo) . " for $gs\n";

  my $h_ref_g = Sets::getIndex("/Users/olivier/PROGRAMS/PLASMODIUM/DATA/GO/GENESETS/$gs.txt");

  open OUT1, ">$ARGV[0].$gs";
  # open matches for each motif

  foreach my $m (@mo) {

    $ta->loadFile("$m.matches");
    $ta->sortbycol(3, 0);

    my $a_ref = $ta->getArray();

    my $cnt = 0;
    foreach my $r (@$a_ref) {
      #print join("\t", @$r) . "\n";
      # if there is a match      
      if (defined($h_ref_g->{$r->[0]})) {
	print OUT1 "$m\t" . join("\t", @$r) . "\n";
      }
      $cnt ++;
      if ($cnt == 1000) {
	last;
      }
    }       

  }
  close OUT1;

  my $genelist = "/Users/olivier/PROGRAMS/PLASMODIUM/DATA/GO/GENESETS/$gs.txt";
  if ($gs =~ /PfEMP1/) {
    $genelist = "/Users/olivier/PROGRAMS/WM-THRESHOLD-FINDER/VAR_RIF/genes_with_ups.txt";
  }

  system("perl ~/PROGRAMS/MOTIFMAPS/draw_motif_map.pl --profiles=$ARGV[0].$gs --genelist=$genelist --seqlen=2000 --leftlabel=\"-2000bp\" --rightlabel=\"0/ATG\"");

  


}
close IN;






#mac131043:TOP250 olivier$ greph PfM combined_matrix | transpose.pl - | sort_column.pl 1 | tail -1 | columns.pl 0 | sed 's/\.top.*//' > PfMC-2TM

# rifin
#for f in `cat rif_enriched_motifs`; do ff=`\ls ../DATA/$f`; perl ~/PERL_SCRIPTS/MyScanACE_output_matches_for_top_N_gene_matches.pl --scanacefile=$ff --n=250 --genelist=../../PLASMODIUM/DATA/GO/GENESETS/rifin.txt --name="genelist_(.+)_pf_u" ; done | tee rif_enriched_motifs_matches.txt

#

       
