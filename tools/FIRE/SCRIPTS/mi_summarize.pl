use lib "$ENV{FIREDIR}/SCRIPTS";


use Table;
use Sets;

use strict;

if (@ARGV == 0) {
  die "Usage: perl mi_summarize.pl -statfile_dna FILE -optimfile_dna FILE -distfile_dna FILE -consfile_dna FILE [ sam for rna] -removecols_draw 0/1 -removecols_exp 0/1 -quantized 0/1 -outsummary FILE -outlabels FILE(colnames) -outsummary FILE -outmatrix FILE -motifnames FILE -outexpfile FILE -outlevel INT\n

-statfile_dna    output from mi_signif
-optimfile_dna   output from mi_optimize (stability)
-distfile_dna    output from mi_dist (position, orientation biases, etc)
-outexpfile      if -removecols is set to 1, output expfile striped from clusters that are not significant 
-expfile_dna     necessary for above
-expfile_rna     same
";
}


my $progdir = "$ENV{FIREDIR}/PROGRAMS";

my $optimfile_dna   = undef;
my $optimfile_rna   = undef;
if (Sets::exist_parameter(\@ARGV, "-optimfile_dna") == 1) {
  $optimfile_dna   = Sets::get_parameter(\@ARGV, "-optimfile_dna");
  if (! -e $optimfile_dna) {
    print "$optimfile_dna does not exist, ignoring.\n";
    $optimfile_dna   = undef; 
  }
}
if (Sets::exist_parameter(\@ARGV, "-optimfile_rna") == 1) {
  $optimfile_rna   = Sets::get_parameter(\@ARGV, "-optimfile_rna");
  if (! -e $optimfile_rna) {
    print "$optimfile_rna does not exist, ignoring.\n";
    $optimfile_rna   = undef; 
  }
}



my $statfile_dna   = undef;
my $statfile_rna   = undef;
if (Sets::exist_parameter(\@ARGV, "-statfile_dna") == 1) {
  $statfile_dna   = Sets::get_parameter(\@ARGV, "-statfile_dna");
  if (! -e $statfile_dna) {
    print "$statfile_dna does not exist, ignoring.\n";
    $statfile_dna   = undef; 
  }
}
if (Sets::exist_parameter(\@ARGV, "-statfile_rna") == 1) {
  $statfile_rna   = Sets::get_parameter(\@ARGV, "-statfile_rna");
  if (! -e $statfile_rna) {
    print "$statfile_rna does not exist, ignoring.\n";
    $statfile_rna   = undef; 
  }
}

my $motifrepfile_dna   = "$statfile_dna.motifs.rep";
my $motifrepfile_rna   = "$statfile_rna.motifs.rep";
if (Sets::exist_parameter(\@ARGV, "-motifrepfile_dna") == 1) {
  $motifrepfile_dna   = Sets::get_parameter(\@ARGV, "-motifrepfile_dna");
}
if (Sets::exist_parameter(\@ARGV, "-motifrepfile_rna") == 1) {
  $motifrepfile_rna   = Sets::get_parameter(\@ARGV, "-motifrepfile_rna");
}

my $distfile_dna   = undef;
my $distfile_rna   = undef;
if (Sets::exist_parameter(\@ARGV, "-distfile_dna") == 1) {
  $distfile_dna   = Sets::get_parameter(\@ARGV, "-distfile_dna");
  if (! -e $distfile_dna) {
    print "$distfile_dna does not exist, ignoring.\n";
    $distfile_dna   = undef; 
  }
}
if (Sets::exist_parameter(\@ARGV, "-distfile_rna") == 1) {
  $distfile_rna   = Sets::get_parameter(\@ARGV, "-distfile_rna");
  if (! -e $distfile_rna) {
    print "$distfile_rna does not exist, ignoring.\n";
    $distfile_rna   = undef; 
  }
}

my $consfile_dna   = undef;
my $consfile_rna   = undef;
if (Sets::exist_parameter(\@ARGV, "-consfile_dna") == 1) {
  $consfile_dna   = Sets::get_parameter(\@ARGV, "-consfile_dna");
  if (! -e $consfile_dna) {
    print "$consfile_dna does not exist, ignoring.\n";
    $consfile_dna   = undef; 
  }
}
if (Sets::exist_parameter(\@ARGV, "-consfile_rna") == 1) {
  $consfile_rna   = Sets::get_parameter(\@ARGV, "-consfile_rna");
  if (! -e $consfile_rna) {
    print "$consfile_rna does not exist, ignoring.\n";
    $consfile_rna   = undef; 
  }
}

my $removecols_draw      = 1;  # for the drawing
if (Sets::exist_parameter(\@ARGV, "-removecols_draw") == 1) {
  $removecols_draw       = Sets::get_parameter(\@ARGV, "-removecols_draw");
}

my $rootdir        = ".";
my $quantized      = Sets::get_parameter(\@ARGV, "-quantized");

my $sortrowsbyphase = 0;
if (Sets::exist_parameter(\@ARGV, "-sortrowsbyphase") == 1) {
  $sortrowsbyphase = Sets::get_parameter(\@ARGV, "-sortrowsbyphase");
}


my $minrobustness = 0;
if (Sets::exist_parameter(\@ARGV, "-minrobustness") == 1) {
  $minrobustness = Sets::get_parameter(\@ARGV, "-minrobustness");
}


my $minp = undef;
if (Sets::exist_parameter(\@ARGV, "-minp") == 1) {
  $minp = Sets::get_parameter(\@ARGV, "-minp");
}


my $minz = undef;
if (Sets::exist_parameter(\@ARGV, "-minz") == 1) {
  $minz = Sets::get_parameter(\@ARGV, "-minz");
}

my $maxrobustness = 0;
if (Sets::exist_parameter(\@ARGV, "-maxrobustness") == 1) {
  $maxrobustness = Sets::get_parameter(\@ARGV, "-maxrobustness");
}

my $outlevel = 1;
if (Sets::exist_parameter(\@ARGV, "-outlevel") == 1) {
  $outlevel = Sets::get_parameter(\@ARGV, "-outlevel");
}

my $outcolumns     = Sets::get_parameter(\@ARGV, "-outcolumns");
if (!defined($outcolumns)) {
  die "Please define the -outcolumns parameter (file in which to store the relevant columns IDs.\n";
}

my $outgoodcolumns = Sets::get_parameter(\@ARGV, "-outgoodcolumns");
if (!defined($outgoodcolumns)) {
  die "Please define the -outgoodcolumns parameter (file in which to store the columns for which there are over-represented motifs.\n";
}
my $outsummary        = Sets::get_parameter(\@ARGV, "-outsummary");
my $outmatrix         = Sets::get_parameter(\@ARGV, "-outmatrix");
my $outfullmatrix     = Sets::get_parameter(\@ARGV, "-outfullmatrix");

my $outdensitymatrix = undef;
if (Sets::exist_parameter(\@ARGV, "-outdensitymatrix") == 1) {
  $outdensitymatrix  = Sets::get_parameter(\@ARGV, "-outdensitymatrix");
}

my $motifnames = undef;
if (Sets::exist_parameter(\@ARGV, "-motifnames") == 1) {
  $motifnames     = Sets::get_parameter(\@ARGV, "-motifnames");
} 

if (Sets::exist_parameter(\@ARGV, "-outnames") == 1) {
  $motifnames     = Sets::get_parameter(\@ARGV, "-outnames");
}

#
#  classical expression vectors
#
my $expfile_dna   = undef;
my $expfile_rna   = undef;
if (Sets::exist_parameter(\@ARGV, "-expfile_dna") == 1) {
  $expfile_dna   = Sets::get_parameter(\@ARGV, "-expfile_dna");
}
if (Sets::exist_parameter(\@ARGV, "-expfile_rna") == 1) {
  $expfile_rna   = Sets::get_parameter(\@ARGV, "-expfile_rna");
}

my $outexpfile     = undef;
if (Sets::exist_parameter(\@ARGV, "-outexpfile") == 1) {
  $outexpfile = Sets::get_parameter(\@ARGV, "-outexpfile");
}


#
#  binary expression matrix
#
my $expfile_b_dna   = undef;
my $expfile_b_rna   = undef;
if (Sets::exist_parameter(\@ARGV, "-expfile_b_dna") == 1) {
  $expfile_b_dna   = Sets::get_parameter(\@ARGV, "-expfile_b_dna");
}
if (Sets::exist_parameter(\@ARGV, "-expfile_b_rna") == 1) {
  $expfile_b_rna   = Sets::get_parameter(\@ARGV, "-expfile_b_rna");
}

my $outexpfile_b     = undef;
if (Sets::exist_parameter(\@ARGV, "-outexpfile_b") == 1) {
  $outexpfile_b = Sets::get_parameter(\@ARGV, "-outexpfile_b");
}


my $removecols_exp      = 0;
if (Sets::exist_parameter(\@ARGV, "-removecols_exp") == 1) {
  $removecols_exp       = Sets::get_parameter(\@ARGV, "-removecols_exp");
}


my $oribiasonly      = undef;
if (Sets::exist_parameter(\@ARGV, "-oribiasonly") == 1) {
  $oribiasonly       = Sets::get_parameter(\@ARGV, "-oribiasonly");
}


if (defined($outexpfile) && !defined($expfile_rna) && !defined($expfile_dna)) {
  die "Please provide an expfile.\n";
}


my $ta = Table->new;

my $h_ref_optim_dna = undef;
if (defined($optimfile_dna)) {
  $ta->loadFile($optimfile_dna);
  $h_ref_optim_dna = $ta->getIndexKV(0, 2);
}
my $h_ref_optim_rna = undef;
if (defined($optimfile_rna)) {
  $ta->loadFile($optimfile_rna);
  $h_ref_optim_rna = $ta->getIndexKV(0, 2);
}

#
#  parse cons file
#

my %cons_dna = ();
if (defined($consfile_dna)) {
  print "Reading DNA conservation file ...";
  $ta->loadFile($consfile_dna);
  my $a_ref_cons = $ta->getArray();
  my $rank = 1;
  my $n = @$a_ref_cons;
  foreach my $r (@$a_ref_cons) {
    $cons_dna{ $r->[0] } = sprintf("%2.2f", $r->[1]/100);
    
    if (defined($r->[2]) && ($r->[2] >= 0.95)) {
      $cons_dna{ $r->[0] } .= "*";
    }
    
    $rank ++;
  }
  print "Done.\n";
}

my %cons_rna = ();
if (defined($consfile_rna)) {
  $ta->loadFile($consfile_rna);
  my $a_ref_cons = $ta->getArray();
  my $rank = 1;
  foreach my $r (@$a_ref_cons) {
    $cons_rna{ $r->[0] } = sprintf("%2.2f", $r->[1]/100);

    if (defined($r->[2]) && ($r->[2] >= 0.95)) {
      $cons_rna{ $r->[0] } .= "*";
    }
    
    $rank ++;
  }
}



#
#  parse distance filess
#
my %dist_dna = ();
my %dist_rna = ();
my %orie_dna = ();
my %orie_rna = ();
my %copy_dna = ();
my %copy_rna = ();

if (defined($distfile_dna)) {
  print "Reading DNA motifs mi_dist output ..."; 
  parse_distance_file($distfile_dna, \%dist_dna, \%orie_dna, \%copy_dna);
  print "Done.\n";
}
if (defined($distfile_rna)) {
  print "Reading RNA motifs mi_dist output ..."; 
  parse_distance_file($distfile_rna, \%dist_rna, \%orie_rna, \%copy_rna);
  print "Done.\n";
}




#
#  parse MOTIF REP files !
#

my %pvalues_dna   = ();
my %pvalues_rna   = ();
my %densities_dna = ();
my %densities_rna = ();
my @labels_dna    = ();
my @labels_rna    = ();
my %ebins_dna     = ();
my %ebins_rna     = ();


if (defined($statfile_dna)) {
  print "Parsing DNA motif file ... ";
  parse_motif_rep_file($motifrepfile_dna, $quantized, \%pvalues_dna, \%densities_dna, \@labels_dna, \%ebins_dna);
  print "Done.\n";
}
if (defined($statfile_rna)) {
  print "Parsing RNA motif file ... ";
  parse_motif_rep_file($motifrepfile_rna, $quantized, \%pvalues_rna, \%densities_rna, \@labels_rna, \%ebins_rna);
  print "Done.\n";
}

if ((@labels_dna == 0) && (@labels_rna > 0)) {
  @labels_dna = @labels_rna;
}


my %max_pv_dna = ();
my %max_pv_rna = ();

if (defined($minp)) {

  if (defined($statfile_dna)) {
    foreach my $k (keys(%pvalues_dna)) {
      my $max = 0;
      foreach my $p (@{ $pvalues_dna{ $k } }) {
	if ($p > $max) {
	  $max = $p;
	}
      }
      $max_pv_dna{ $k } = $max;
      #print "$k\t$max\n";
    }
  }
  
  if (defined($statfile_rna)) {
    foreach my $k (keys(%pvalues_rna)) {
      my $max = 0;
      foreach my $p (@{ $pvalues_rna{ $k } }) {
	if ($p > $max) {
	  $max = $p;
	}
      }
      $max_pv_rna{ $k } = $max;
      #print "$k\t$max\n";
    }
  }
  
  
}


#
#  parse stat file
#
my @MOT = ();

my $dec = 0;
if ($quantized == 0) {
  $dec = 1;
}

my @SIGSETS = ("OK");
if ($outlevel == 2) {
  push @SIGSETS, "OK-NO-SEED";
}

if (defined($statfile_dna)) {
  print "Reading mi_signif output file for DNA motifs ...";

  $ta->loadFile($statfile_dna);
  my $a_ref = $ta->getArray();
  foreach my $r (@$a_ref) {
    
    #print "Test $r->[0]\t$r->[4]\t" . $max_pv_dna{$r->[0]} . "\n";
    
    if (Sets::in_array($r->[5], @SIGSETS) && ($r->[4] >= $minrobustness))  {

      # added filter by p-value
      next if (defined($minp) && ($max_pv_dna{$r->[0]} < $minp));

      next if (defined($minz) && ($r->[3] < $minz));

      
      my %a = ("MOTIF" => $r->[0], 
	       "COPIES"=> $copy_dna{ $r->[0] }, 
	       "MI"    => $r->[1], 
	       "RNA"   => 0, 
	       "Z"     => $r->[3], 
	       "S"     => $h_ref_optim_dna->{$r->[0]}, 
	       "RANK"  => $r->[2], 
	       "R"     => $r->[4],
	       "SEED"  => $r->[6+$dec], 
	       "DIST"  => $dist_dna{ $r->[0] }, 
	       "ORIE"  => $orie_dna{ $r->[0] }, 
	       "CONS"  => (defined($consfile_dna)?$cons_dna{ $r->[0] }:"-"));
      push @MOT, \%a;

    }
  }
  print "Done.\n";
  
  print "Got " . scalar(@MOT) . " motifs\n";
}


if (defined($statfile_rna)) {

  print "Reading mi_signif output file for RNA motifs ...";

  $ta->loadFile($statfile_rna);
  my $a_ref = $ta->getArray();
  
  foreach my $r (@$a_ref) {
    if (Sets::in_array($r->[5], @SIGSETS) && ($r->[4] >= $minrobustness)) {

      # added filter by p-value
      next if (defined($minp) && ($max_pv_rna{$r->[0]} < $minp));

      next if (defined($minz) && ($r->[3] < $minz));

      next if (defined($oribiasonly) && ($oribiasonly == 1) && ($orie_rna{ $r->[0] } == 0));

      my %a = ("MOTIF" => $r->[0], 
	       "COPIES"=> $copy_rna{ $r->[0] }, 
	       "MI"    => $r->[1], 
	       "RNA"   => 1, 
	       "Z"     => $r->[3], 
	       "S"     => $h_ref_optim_rna->{$r->[0]}, 
	       "RANK"  => $r->[2], 
	       "R"     => $r->[4], 
	       "SEED"  => $r->[6+$dec], 
	       "DIST"  => $dist_rna{ $r->[0] }, 
	       "ORIE"  => $orie_rna{ $r->[0] }, 
	       "CONS"  => (defined($consfile_rna)?$cons_rna{ $r->[0] }:"-"));

      push @MOT, \%a;
    }
  }
  print "Done.\n";

  print "Got " . scalar(@MOT) . " motifs\n";

}



#
#  sort motifs by MI (default), then build A_MOT
#
@MOT = sort { $b->{MI} <=> $a->{MI} } @MOT;
my @A_MOT = ();
foreach my $r (@MOT) {
  push @A_MOT, $r->{MOTIF};
}




#
# remove pvalues for non-signif motifs, merge both matrices at the same time 
#
my %pvalues = ();
foreach my $k (keys(%pvalues_dna)) {
  if (Sets::in_array($k, @A_MOT)) {
    $pvalues{ $k } = $pvalues_dna{ $k };
  }
}

foreach my $k (keys(%pvalues_rna)) {
  if (Sets::in_array($k, @A_MOT)) {
    $pvalues{ $k } = $pvalues_rna{ $k };
  }
}

#
# same for densities
#
my %densities = ();
foreach my $k (keys(%densities_dna)) {
  if (Sets::in_array($k, @A_MOT)) {
    $densities{ $k } = $densities_dna{ $k };
  }
}

foreach my $k (keys(%densities_rna)) {
  if (Sets::in_array($k, @A_MOT)) {
    $densities{ $k } = $densities_rna{ $k };
  }
}

#
# remove BAD columns
#

# 1) first build an index of the good columns
my %GOODCLU = ();
foreach my $k (keys(%ebins_dna)) {
  next if (!Sets::in_array($k, @A_MOT));
  my $r = $ebins_dna{$k};
  foreach my $s (@$r) {
    $GOODCLU{$s} = 1;
  }
}
foreach my $k (keys(%ebins_rna)) {
  next if (!Sets::in_array($k, @A_MOT));
  my $r = $ebins_rna{$k};
  foreach my $s (@$r) {
    $GOODCLU{$s} = 1;
  }
}


my $n = @labels_dna;

my @columnstoshow = ();
for (my $j=0; $j<$n; $j++) {
  push @columnstoshow, $j;
}


my @goodcolumns = ();
for (my $j=0; $j<$n; $j++) {
  if (defined($GOODCLU{$j})) {
    push @goodcolumns, $j;
  }
}


#
# output columns IDs
#
print "Writing good columns to $outgoodcolumns ...";
Sets::writeSet(\@goodcolumns, $outgoodcolumns);
print "Done.\n";


#
# output FULL p-value matrix
#
open OUT, ">$outfullmatrix";
print OUT "\t" . join("\t", @labels_dna) . "\n";
foreach my $k (@A_MOT) {
  print OUT "$k\t"; print OUT join("\t", @{ $pvalues{$k} }); print OUT "\n";
}
close OUT;



if ($removecols_draw == 1) {


  # 2) remove bad columns from the labels
  my @a           = ();
  @columnstoshow  = (); # reset, we need to reduce it
  for (my $j=0; $j<$n; $j++) {
    if (defined($GOODCLU{$j})) {
      push @columnstoshow, $j;
      push @a, $labels_dna[$j];
    }
  }
  @labels_dna = @a;
  
  
  # 3) remove bad columsn from the pvalue matrix
  
  foreach my $k (@A_MOT) {
    my $r = $pvalues{$k};
    
    my @a = ();
    for (my $j=0; $j<@$r; $j++) {
      if (defined($GOODCLU{$j})) {
	push @a, $r->[$j];
      }
    }
    $pvalues{$k} = \@a;

  }

  # same for density matrix
  foreach my $k (@A_MOT) {
    my $r = $densities{$k};
    
    my @a = ();
    for (my $j=0; $j<@$r; $j++) {
      if (defined($GOODCLU{$j})) {
	push @a, $r->[$j];
      }
    }
    $densities{$k} = \@a;

  }
}


my %GENES_INEXP = ();
if (defined($outexpfile)) {
  
  print "Creating new expression file with only good clusters/bins ...";
  
  my $he = undef;
  my @INEXP = ();
  
  if ((defined($expfile_dna) && !defined($expfile_rna)) || (defined($expfile_rna) && !defined($expfile_dna))) {
    
    my $expfile = undef;
    if (defined($expfile_dna)) {
      $expfile = $expfile_dna;
    }
    if (defined($expfile_rna)) {
      $expfile = $expfile_rna; 
    }
    
    $ta->loadFile($expfile);
    my $a_ref = $ta->getArray();
    $he = shift @$a_ref;
    
    foreach my $r (@$a_ref) {
      push @INEXP, $r;  # do it, no filter
    }
    
  } elsif (defined($expfile_rna) && defined($expfile_dna)) {
    
    $ta->loadFile($expfile_dna);
    my $a_ref_dna = $ta->getArray();
    $ta->loadFile($expfile_rna);
    my $h_ref_rna = $ta->getIndexKV(0,1);
    $he = shift @$a_ref_dna;
    
    
    foreach my $r (@$a_ref_dna) {
      if (defined($h_ref_rna->{$r->[0]})) {
	push @INEXP, $r;
      }
    }
    
  } else {
    die "WTF?\n";
  }
  
  
  #
  #  reset expression bin counters if quantized
  #
  
  open OUT, ">$outexpfile" or die "Cannot open outexpfile: $outexpfile.\n";
  print OUT join("\t", @$he) . "\n";
  foreach my $r (@INEXP) {
      print OUT join("\t", @$r) . "\n";
      $GENES_INEXP{ $r->[0] } = 1;
  }
  close OUT;
  print "Done.\n";
  
}


#
# order rows by phase
#
if ($sortrowsbyphase == 1) {
  
  print "Sorting motifs by phase ... ";

  # determine the position of the peak
  my @PE = ();
  foreach my $m (@A_MOT) {
    my $max_p = -100000;
    my $max_i =  undef;
    for (my $i=0; $i<@{ $pvalues{$m} }; $i++) {
      if ($pvalues{$m}->[$i] > $max_p) {
	$max_p = $pvalues{$m}->[$i];
	$max_i = $i;
      }
    }
    my @a = ($m, $max_i);
    push @PE, \@a;
  }
  
  @PE = sort { $a->[1] <=> $b->[1] } @PE;

  @A_MOT = ();
  foreach my $r (@PE) {
    print "$r->[0]\t$r->[1]\n";
    push @A_MOT, $r->[0];
  }
  
  print "Done.\n";
  
}



#
#  START:  output SUMMARY file
#

print "Outputing summary ... ";
open OUT, ">$outsummary" or die "cannot open $outsummary\n";

foreach my $m (@A_MOT) {

  my $thei = undef;
  for (my $i=0; $i<@MOT; $i++) {
    if ($MOT[$i]->{MOTIF} eq $m) {
      $thei = $i; last;
    }
  }
  
  my $r = $MOT[$thei];
  
  print OUT "$r->{MOTIF}\t$r->{RNA}\t$r->{COPIES}\t$r->{MI}\t$r->{RANK}\t$r->{Z}\t$r->{R}/$maxrobustness\t$r->{S}\t$r->{SEED}\t$r->{DIST}\t$r->{ORIE}\t$r->{CONS}";

  #
  # add enriched clusters
  #
  if (defined($ebins_dna{ $r->{MOTIF} }) && (@{ $ebins_dna{ $r->{MOTIF} } } > 0)) {
    print OUT "\t" . join("\t", @{ $ebins_dna{ $r->{MOTIF} } });
  } elsif (defined($ebins_rna{ $r->{MOTIF} }) && (@{ $ebins_rna{ $r->{MOTIF} } } > 0)) {
    print OUT "\t" . join("\t", @{ $ebins_rna{ $r->{MOTIF} } });
  }

  print OUT "\n";
  
}
close OUT;
print "Done.\n";

#
# STOP: output SUMMARY file
#


#
# START: input FULL binary file, then output BINARY with columns in same order as 
#

# 1. input full DNA binary
if (defined($expfile_b_dna) || defined($expfile_b_rna)) {
  
  #die "Please define -outexpfile.\n" if (!defined($outexpfile));

  my %HE = ();  # will contain all expression
  
  if (defined($expfile_b_dna)) {
    $ta->loadFile($expfile_b_dna);
    my $a_ref_expb = $ta->getArray();  
    my $r_m = shift @$a_ref_expb;  # = motif list
    shift @$r_m;
    foreach my $r (@$a_ref_expb) {
      my $g = shift @$r;
      next if (defined($outexpfile) && !defined($GENES_INEXP{$g}));
      for (my $i=0; $i<@$r; $i++) {
	die "bug?\n" if (!defined($r_m->[$i]));
	$HE{ $g } { $r_m->[$i] } = $r->[$i];
      }	
    }
  }
  
  if (defined($expfile_b_rna)) {
    $ta->loadFile($expfile_b_rna);
    my $a_ref_expb = $ta->getArray();  
    my $r_m = shift @$a_ref_expb;  # = motif list
    shift @$r_m;
    foreach my $r (@$a_ref_expb) {
      my $g = shift @$r;
      next if (defined($outexpfile) && !defined($GENES_INEXP{$g}));
      for (my $i=0; $i<@$r; $i++) {
	die "bug?\n" if (!defined($r_m->[$i]));
	$HE{ $g } { $r_m->[$i] } = $r->[$i];
      }	
    }
  }

  
  die "please specify -outexpfile_b.\n" if (!defined($outexpfile_b));
  open OUTG, ">$outexpfile_b" or die "Cannot open outexpfile_b $outexpfile_b for writing.\n";
  

  foreach my $m (@A_MOT) {
    print OUTG "\t$m";
  }
  print OUTG "\n";

  foreach my $g (keys(%HE)) {
    my $h = $HE{$g};
    print OUTG "$g";
    foreach my $m (@A_MOT) {
      print OUTG "\t" . $h->{$m};
    }
    print OUTG "\n";
  }
  
  close OUTG;

  print "$outexpfile_b created.\n";

}


# 2. input full RNA binary



# 3. output full DNA/RNA binary




#
# output columns IDs to show
#
Sets::writeSet(\@columnstoshow, $outcolumns);

  

#
# output motif names
#
my @A_MOT_BIS = @A_MOT;
foreach my $r (@A_MOT_BIS) {
  $r .= "\t-";
}
Sets::writeSet(\@A_MOT_BIS, $motifnames);



#
# output matrix
#
open OUT, ">$outmatrix";
print OUT "\t" . join("\t", @labels_dna) . "\n";
foreach my $k (@A_MOT) {
  print OUT "$k\t"; print OUT join("\t", @{ $pvalues{$k} }); print OUT "\n";
}
close OUT;


if (defined($outdensitymatrix)) {
  open OUT, ">$outdensitymatrix";
  print OUT "\t" . join("\t", @labels_dna) . "\n";
  foreach my $k (@A_MOT) {
    print OUT "$k\t"; print OUT join("\t", @{ $densities{$k} }); print OUT "\n";
  }
  close OUT;
}


sub parse_motif_rep_file {
  
  my ($repfile, $quantized, $h_ref_pvalues, $h_ref_densities, $a_ref_labels, $h_ref_ebins) = @_;



  my @MOTIFS      = ();
  my $cnt         = 1;
  my %FREQ        = ();
  my @N           = ();
  my %avg_density = ();
  my %sample_size = ();

  #
  # read in motif repfile
  #
  
  my $ta = Table->new;
  
  $ta->loadFile($repfile);
  my $a_ref = $ta->getArray();
  
  foreach my $r (@$a_ref) {
    
    #
    #  motif is always first column
    #
    push @MOTIFS, $r->[0] if !Sets::in_array($r->[0], @MOTIFS);
    
    if ($quantized == 1) {
      push @{ $FREQ{ $r->[0] } }, $r->[4];
      push @N, $r->[5];  # number of item per cluster
      
      $sample_size{ $r->[0] } += $r->[5];
      $avg_density{ $r->[0] } += int( 0.5 + $r->[4] * $r->[5] );
      
    } else {
      push @{ $FREQ{ $r->[0] } }, $r->[6];
      push @N, $r->[7];  # number of item per cluster
      
      $sample_size{ $r->[0] } += $r->[7];
      $avg_density{ $r->[0] } += int( 0.5 + $r->[6] * $r->[7] );
    }
    
    my $label = undef;
    if ($quantized == 1) {
      $label = "C$r->[1]";
    } else {
      $label = "[$r->[4];$r->[5]]";
    }
    
    push @$a_ref_labels, $label if (!Sets::in_array($label, @$a_ref_labels));

    $cnt ++;
  }
  

  #
  # calculate average density
  #
  foreach my $k (keys(%avg_density)) {
    $avg_density{$k} = $avg_density{$k} / $sample_size{$k};
  }
  
  
  
  #
  # record significant exp bins for each motif
  #

  
  my $cnt = 0;
  foreach my $k (@MOTIFS) {
    
    #print "Processing $k.                        \n";
    
    my $n = scalar( @{ $FREQ{$k} } ); 
    
    #print "I have $n bins\n";

    for (my $i=0; $i<$n; $i++) {
      
      

      my $bk = undef;
      my $bn = undef;
      
      my $the_ebin = undef;
      
      $bk = int( 0.5 + $FREQ{$k}->[$i] * $N[$i] );
      $bn = $N[$i];
      $the_ebin = $i;
      
      my $bN = $avg_density{$k};    
      my $p2 = undef;
      
      my $p1 = binom_test_greater( $bk, $bn, $bN, \$p2);
      
      my $p  = undef;
      
      if ($p1 < 1e-100) {
	$p1 = 1e-100;
      }
      if ($p2 < 1e-100) {
	$p2 = 1e-100;
      }
      if ($p1 < $p2) {
	$p = sprintf("%4.3f", - Sets::log10($p1));
      } else {
	$p = sprintf("%4.3f", + Sets::log10($p2));
      }
      
      if ($p1 < 0.05 / $n) {
	push @{ $h_ref_ebins->{ $k } }, $the_ebin;
      }
      
      $h_ref_pvalues   ->{ $k }->[ $i ] = $p;

      $h_ref_densities ->{ $k }->[ $i ] = $FREQ{$k}->[$i];
      
    }
    
    
    $cnt ++;    
 
    #print "Done\n";
  }
  
}



sub binom_test_greater {
  
  my ($n, $k, $p, $p1) = @_;

  my $todo = "$progdir/binom_test_greater $n $k $p 1";
  my $out = `$todo`;  
  $out =~ s/[\n\r]//g;

  my @a = split /\t/, $out, -1;
  
  $out = $a[0];

  $$p1 = $a[1];
  
  return $out;

}



sub parse_distance_file {
  my ($distfile, $RDIST, $RORIE, $RCOPY) = @_;
  
  my $o_t = 1;

  #
  #  load distance info
  #
  $ta->loadFile($distfile);
  my $a_ref_dist = $ta->getArray();

  my %TMPDIST   = ();
  my @MOTIFS    = ();
  foreach my $r (@$a_ref_dist) {
    my $mo = shift @$r; 
    push @MOTIFS, $mo if (!Sets::in_array($mo, @MOTIFS));
    my $ke = shift @$r;
    $TMPDIST{$mo}{$ke} = $r;
  }
  
  foreach my $m (@MOTIFS) {

    if ($TMPDIST{ $m }{"d_avg"}->[2] <= 100) {    
      $RDIST->{ $m } = 1;
    } else {
      $RDIST->{ $m } = 0;
    }
    
    if (($TMPDIST{ $m }{"o5"}->[2] == 0 ) && ($TMPDIST{ $m }{"o3"}->[2] >= $o_t)) {
      $RORIE->{ $m } = 1;
    } elsif (($TMPDIST{ $m }{"o5"}->[2] >= $o_t) && ($TMPDIST{ $m }{"o3"}->[2] ==  0)) {
      $RORIE->{ $m } = 2;
    } else {
      $RORIE->{ $m } = 0;
    }
    
    my $i=1; my $bi = undef; my $zz = -10000; my $mm = -10000;
    while (defined($TMPDIST{ $m }{"c$i"})) {
      if (($TMPDIST{ $m }{"c$i"}->[3] > $zz) && ($TMPDIST{ $m }{"c$i"}->[1] > $mm))  {
	$bi = $i; $mm = $TMPDIST{ $m }{"c$i"}->[1]; $zz = $TMPDIST{ $m }{"c$i"}->[3];
      }
      $i++;
    }
    
    $RCOPY->{$m} = $bi;
  }
  
}
