use lib "$ENV{FIREDIR}/SCRIPTS";

use Sets;
use Table;
use Getopt::Long;

use strict;

if (@ARGV == 0) {
  die "Usage: perl combine_multiple_fire_runs.pl --expfiles=FILE[S] --species=STR --exptype=STR [ --tc=FLT --showall=0/1 --twocats=0/1 --dophase1=0/1 --dophase2=0/1 --restrict=FILE  ]\n";
}


my $expfiles = undef;
my $exptype  = undef;
my $rna      = undef;
my $tc       = 0.7;
my $showall  = 0;
my $dophase1 = 1;
my $dophase2 = 0;
my $sortcondsbystage = 0;
my $twocats  = 0;
my $species  = undef;
my $clust    = undef;
my $sortmotifsbyphase = 0;
my $restrict = undef;
my $jn_t_col = 6;   # threshold for including a column in the figure
my $sortcondsbyseries = undef;
my $restrictfile = undef;
my $docons   = 1;
my $mincntrow=0;
my $onlyclust=0;
my $walltime = undef;
my $listexpfiles = undef;
my $verbose      = 0;
my $doclust      = 1;
my $submit       = 1;
my $platform     = undef;

GetOptions ('expfiles=s'             => \$expfiles,
	    'exptype=s'              => \$exptype,
	    'listexpfiles=s'         => \$listexpfiles,
	    'doclust=s'              => \$doclust,
	    'submit=s'               => \$submit,
	    'tc=s'                   => \$tc,
	    'verbose=s'              => \$verbose,
	    'docons=s'               => \$docons,
	    'restrict=s'             => \$restrict,
	    'sortmotifsbyphase=s'    => \$sortmotifsbyphase,
	    'dophase1=s'             => \$dophase1,
	    'dophase2=s'             => \$dophase2,
	    'sortcondsbystage=s'     => \$sortcondsbystage,
	    'twocats=s'              => \$twocats,
	    'species=s'              => \$species,
	    'restrictfile=s'         => \$restrictfile,
	    'clust=s'                => \$clust,
	    'jn_t_col=s'             => \$jn_t_col,
	    'mincntrow=s'            => \$mincntrow,
	    'onlyclust=s'            => \$onlyclust,
	    'walltime=s'             => \$walltime,
	    'platform=s'             => \$platform,
	    'sortcondsbyseries=s'    => \$sortcondsbyseries,
	    'showall=s'              => \$showall);


if (!defined($species)) {
  die "Please specify a species.\n";
}

#if (!defined($exptype)) {
#  die "Please specify exp type (continuous or discrete).\n";
#}

if ($dophase2 == 1) {
  $dophase1 = 0;
}

my $dirname = undef;
if (defined($expfiles)) {
  $dirname = Sets::dirname($expfiles);
} else {
  $dirname = ".";
}

my $h_ref_good = undef;
if (defined($restrict)) {
  $h_ref_good = Sets::getIndex($restrict);
}

if ($dophase1 == 1) {

  #1 clustering

  mkdir "$dirname/META" if (! -e "$dirname/META");

  if ($doclust == 1) {
    
    #
    # cluster DNA motifs
    #
    print "Clustering DNA motifs.\n";
    
    my $todo = "perl $ENV{FIREDIR}/MORESCRIPTS/cluster_motifs_from_multiple_fire_runs.pl ";
    if (defined($expfiles)) {
      $todo .= " --expfiles=\"$expfiles\" ";
    } elsif (defined($listexpfiles)) {
      $todo .= " --listexpfiles=$listexpfiles ";
    }
    $todo .= " --tc=$tc --rna=0 --outfile=$dirname/META/meta-motifs_dna --outfull=$dirname/META/meta-motifs_dna_full ";
    
    if ($verbose == 1) {
      print "$todo\n";
    }
    system($todo);

    #
    # cluster RNA motifs
    #
    print "Clustering RNA motifs.\n";
    $todo  = "perl $ENV{FIREDIR}/MORESCRIPTS/cluster_motifs_from_multiple_fire_runs.pl ";
    if (defined($expfiles)) {
      $todo .= " --expfiles=\"$expfiles\" ";
    } elsif (defined($listexpfiles)) {
      $todo .= " --listexpfiles=$listexpfiles ";
    }
    $todo .= " --tc=$tc --rna=1 --outfile=$dirname/META/meta-motifs_rna --outfull=$dirname/META/meta-motifs_rna_full ";
    if ($verbose == 1) {
      print "$todo\n";
    }
    
    system($todo);
        

  }

  system("cat $dirname/META/meta-motifs_dna");
  print "\n";
  system("cat $dirname/META/meta-motifs_rna");
  print "Results are OK ? press enter to continue, ctrl-C to stop.\n";
  <STDIN>;

  #2 re-run FIRE on motif clusters

  if (defined($expfiles)) {

    print "Copying expfiles to META/ directory.\n";
    my $todo = "cp $expfiles $dirname/META/";
    system($todo);

    if ($onlyclust == 0) {
      print "Running FIRE in non-discovery mode on the cluster representatives.\n";
      $todo = "perl $ENV{FIREDIR}/fire.pl --expfiles=\"$dirname/META/*.*\" --exptype=$exptype --species=$species --doskipdiscovery=1 --docons=0 --dogomotifs=0 --dogoclusters=0 --domicombine=0 --motiffile_dna=$dirname/META/meta-motifs_dna --motiffile_rna=$dirname/META/meta-motifs_rna --jn_t=0 --submit=$submit ";
      if (defined($platform)) {
	$todo .= " --platform=$platform ";
      }
      if (defined($walltime)) {
	$todo .= " --walltime=$walltime ";
      }

      $todo .= " > /dev/null";

      if ($verbose == 1) {
	print "$todo\n";
      }

      system($todo);

    }

  } else {
    
    my $ta     = Table->new;
    $ta->loadFile($listexpfiles);
    my $a_ref_f = $ta->getArray();

    foreach my $r (@$a_ref_f) {
      system("cp $r->[0] $dirname/META/");
      
      my $ff = Sets::filename($r->[0]);

      my $todo = "perl $ENV{FIREDIR}/fire.pl --expfile=$dirname/META/$ff --exptype=$r->[1] --species=$species --doskipdiscovery=1 --docons=0 --dogomotifs=0 --dogoclusters=0 --domicombine=0 --motiffile_dna=$dirname/META/meta-motifs_dna --motiffile_rna=$dirname/META/meta-motifs_rna --jn_t=0 --submit=$submit ";
      if (defined($platform)) {
	$todo .= " --platform=$platform ";
      }
      if (defined($walltime)) {
	$todo .= " --walltime=$walltime ";
      }

      if ($verbose == 0) {
	$todo .= " > /dev/null";
      }
      system($todo);
    }
    
  }
}

if ($dophase2 == 1) {

  #3 gather info from .summary files
  my $ta     = Table->new;
  my $a_ref_f = Sets::getFiles("$dirname/META/*.txt $dirname/META/*.kgg $dirname/META/*.logratios");
  
  my %DNA = ();
  my %RNA = ();
  my %CON = ();
  
  #
  # 
  #
  print "Read " . scalar(@$a_ref_f) . " files\n";

  my $h_ref_exptypes = {};

  if (!defined($exptype)) {
    $ta->loadFile($listexpfiles);
    my $a_ref_l = $ta->getArray();
    foreach my $r (@$a_ref_l) {
      my $m = Sets::filename($r->[0]);
      $h_ref_exptypes->{ "$dirname/META/$m" } = $r->[1]; 
    }
  }

  my %NAMES      = ();
  my %CONSINDEX  = ();
  foreach my $f (@$a_ref_f) {

    if ($verbose == 1) {
      print "ANALYZING $f.\n";
    }

    my $fi  = Sets::filename($f);
    my $ff  = "$f\_FIRE/DNA_RNA/$fi.summary";

    my $fo1 = "$f\_FIRE/DNA_RNA/$fi.signif.dna";
    my $fo2 = "$f\_FIRE/DNA_RNA/$fi.signif.rna";
    my $fn  = "$f\_FIRE/DNA_RNA/$fi.motifnames";

    if (-e $ff) {

      my $h_ref_dir = undef;

      
      if (!defined($exptype) && !defined($h_ref_exptypes->{$f})) {
	die "problem .. $f ?\n";
      }

      if ((defined($exptype) && ($exptype eq "continuous")) || (!defined($exptype) && ($h_ref_exptypes->{$f} eq "continuous"))) {

	print "Entered here .. $f \n";
	$ta->loadFile($fo1);
	my $a_ref_dir = $ta->getArray();
	foreach my $o (@$a_ref_dir) {
	  $h_ref_dir->{$o->[0]} = $o->[6];
	}
	$ta->loadFile($fo2);
	$a_ref_dir = $ta->getArray();
	foreach my $o (@$a_ref_dir) {
	  $h_ref_dir->{$o->[0]} = $o->[6];
	}	
      }

      #	
      # load file name
      #
      $ta->loadFile($fn);
      my $a_ref_names = $ta->getArray();
      foreach my $r (@$a_ref_names) {
	$NAMES{ $r->[0] } = $r->[1];
	#print "$r->[0]\t$r->[1]\t$ff\n";
      }

      $ta->loadFile($ff);
      my $a_ref = $ta->getArray();
    
      #
      # go over all motifs
      #
      foreach my $r (@$a_ref) {	

	print "$r->[0]\t$r->[1]\n";

	# skip motif in twocats case if not enriched in C1
	if (($twocats == 1) && (!defined($r->[12]) || ($r->[12] eq "") || ($r->[12] == 0))) {
	  next;
	}
	
	# skip if not in short list
	if (defined($restrict) && !defined($h_ref_good->{$r->[0]})) {
	  next;
	}
	
	

	$CONSINDEX{$r->[0]} = $r->[11];

	my $ro = $r->[6]; $ro =~ s/\/10$//;
	
	my $di = 1;
	if ((defined($exptype) && ($exptype eq "continuous")) || (!defined($exptype) && ($h_ref_exptypes->{$f} eq "continuous"))) {
	  $di = $h_ref_dir->{ $r->[0] };
	  print "direction = $di ... for $f\n"

	} 

	if ($r->[1] == 0) {
	  $DNA{ $r->[0] }{ $fi } = Sets::sign($di) * $ro;
	} else {
	  $RNA{ $r->[0] }{ $fi } = Sets::sign($di) * $ro;
	}
	
	$CON{$fi} = 1 if ($ro >= $jn_t_col);  # do not include the bad cols
      }
      
      
    } else {
      print "$ff does not exist.\n";
    }
  }

  

  #4 output matrix
  my @MOTIFS = (keys(%DNA), keys(%RNA));

  print "Got " . scalar(@MOTIFS) . " motifs.\n";

  open OUT, ">$dirname/META/matrix-motif-type";
  foreach my $m (keys(%DNA)) {
    print OUT "$m\t0\n";
  }
  foreach my $m (keys(%RNA)) {
    print OUT "$m\t1\n";
  }
  close OUT;


  if ($docons == 1) {
    open OUT, ">$dirname/META/motif-cons";
    foreach my $k (keys(%CONSINDEX)) {
      print OUT "$k\t$CONSINDEX{$k}\n";
    }
    close OUT;
  }



  open OUT, ">$dirname/META/matrix";

  
  my @CONDS1  = ();  # first line, the one that is displayed
  my @CONDS2  = ();  # second line, which contain the file names

  if ($sortcondsbystage == 1) {
    @CONDS1 = sort sortByStage keys(%CON);
    @CONDS2 = @CONDS1;

  } elsif (defined($sortcondsbyseries)){

    # note that in that case, we specifically want all conditions
    my $a_ref_series = undef;
    if (defined($sortcondsbyseries)) {
      $ta->loadFile($sortcondsbyseries);
      $a_ref_series = $ta->getArray();
      foreach my $c (@$a_ref_series) {
	push @CONDS1, $c->[1];
	push @CONDS2, $c->[0];
      }
    }
    
  } elsif (defined($listexpfiles)) {
    my $ta     = Table->new;
    $ta->loadFile($listexpfiles);
    my $a_ref_f = $ta->getArray();
    foreach my $f (@$a_ref_f) {
      push @CONDS1, Sets::filename($f->[0]);
    }
    @CONDS2 = @CONDS1;

  } else {
    @CONDS1 = sort keys(%CON);
    @CONDS2 = @CONDS1;
  }


  foreach my $c (@CONDS1) {
    print OUT "\t$c";
  }
  print OUT "\n";
  foreach my $c (@CONDS2) {
    print OUT "\t$c";
  }
  print OUT "\n";

  #
  # determine if there are some motifs to leave out
  #
  my %DNA_OUT = ();
  foreach my $m (keys(%DNA)) {
    my $cntdna = 0;
    foreach my $c (@CONDS2) {
      $cntdna ++ if ($DNA{$m}{$c} > 0);
    }
    if ($cntdna < $mincntrow) {
      $DNA_OUT{$m} = 1;
    }
  }
  my %RNA_OUT = ();
  foreach my $m (keys(%RNA)) {
    my $cntrna = 0;
    foreach my $c (@CONDS2) {
      $cntrna ++ if ($RNA{$m}{$c} > 0);
    }
    if ($cntrna < $mincntrow) {
      $RNA_OUT{$m} = 1;
    }
  }

  
  #
  # output matrix (only good lines)
  #
  foreach my $m (keys(%DNA)) {
    next if (defined($DNA_OUT{$m}));

    print OUT "$m\t" . $NAMES{$m}; 
    foreach my $c (@CONDS2) {
      print OUT "\t" . (defined($DNA{$m}{$c})?$DNA{$m}{$c}:0);
    }
    print OUT "\n";
  }

  foreach my $m (keys(%RNA)) {
    next if (defined($RNA_OUT{$m}));

    print OUT "$m\t" . $NAMES{$m};
    foreach my $c (@CONDS2) {
      print OUT "\t" . (defined($RNA{$m}{$c})?$RNA{$m}{$c}:0);
    }
    print OUT "\n";
  }
  close OUT;



  #5 draw matrix
  my $todo = "perl $ENV{FIREDIR}/MORESCRIPTS/draw_mi_matrix.pl -matrixfile $dirname/META/matrix -colmap $ENV{FIREDIR}/SCRIPTS/HEATMAPS/cmap_bones.txt -w 10 -motiftypefile $dirname/META/matrix-motif-type ";
  if ($sortmotifsbyphase == 1) {
    $todo .= " -sortmotifsbyphase 1 ";
  } else  {
    
    if (!defined($clust)) {
      $clust = int( 0.5 + sqrt( scalar(@MOTIFS) ));
    }
    
    $todo .= " -clust $clust ";
  }	
  
  if (defined($restrictfile)) {
    $todo .= " -restrictfile $restrictfile ";
  }

  if ($docons == 1) {
    $todo .= " -consfile $dirname/META/motif-cons ";
  }

  if ($twocats == 1) {
    $todo .= " -posonly 1 ";
  }

  print "$todo\n";
  system($todo);


}



sub sortByStage {
  #my ($a, $b) = @_;
  my ($aa) = $a =~ /^stage(\d+)\-/;
  my ($bb) = $b =~ /^stage(\d+)\-/;
  return $aa <=> $bb;
}
