#!/usr/bin/perl

BEGIN{ 

  if ((!$ENV{FIREDIR}) || ($ENV{FIREDIR} eq '')) {
    print "The FIREDIR environment variable is not set. Please set it using export or setenv (see FIRE tutorial online).\n";
    exit;
  }

}

use lib "$ENV{FIREDIR}/SCRIPTS";

use Fire;
use Sets;
use PBS;
use Table;
use Getopt::Long;
use strict;

my $progdir = "$ENV{FIREDIR}/PROGRAMS";
my $pwd     = Sets::getpwd();

if (@ARGV == 0) {
  die "Usage: perl optimize_each_motif_to_WM_using_grid.pl --expfile= --fastafile_dna=FILE --fastafile_rna=FILE\n";
}


my $expfile       = undef;
my $add           = 3;
my $species       = undef;
my $exptype       = undef;
my $fastafile_dna = undef;
my $fastafile_rna = undef;
my $verbose       = 0;
my $walltime      = "00:30:00";
my $submit        = 0;

GetOptions ('expfile=s'       => \$expfile,
	    'exptype=s'       => \$exptype,
	    'species=s'       => \$species,
	    'add=s'           => \$add,
	    'submit=s'        => \$submit,
	    'verbose=s'       => \$verbose,
	    'walltime=s'      => \$walltime,
	    'fastafile_dna=s' => \$fastafile_dna,
	    'fastafile_rna=s' => \$fastafile_rna);


#
# load in .summary file (DNA+RNA)
#

my $expfile_only = Sets::filename($expfile); 
my $summaryfile  = "$expfile\_FIRE/DNA_RNA/$expfile_only.summary";

my $targetdir    = "$summaryfile.optimization_OUT";
mkdir $targetdir if (! -e $targetdir);

if (-e $summaryfile) {
  print "Correctly found $summaryfile.\n";
} else {
  die "Did not find $summaryfile, exiting.\n";
}

my $a_ref_m      = Fire::loadFireMotifSummaryArray($summaryfile);


my $motif_num = 0;

foreach my $r (@$a_ref_m) {

  
  my $fastafile = undef;
  if ($r->{RNA} == 0) {
    
    if (defined($fastafile_dna) && (-e $fastafile_dna)) {
      $fastafile = $fastafile_dna;
    } else {
      die "Cannot find $fastafile_dna.\n";
    }
  } else {
    
    if (defined($fastafile_rna) && (-e $fastafile_rna)) {
      $fastafile = $fastafile_rna;
    } else {
      die "Cannot find $fastafile_rna.\n";
    }
    
  }
  
  my $todo = "$progdir/mi_optimize_motif_WM -motif $r->{MOTIF} -fastafile $fastafile -quantized $exptype -expfile $expfile -verbose $verbose -add $add -rna $r->{RNA} -outfile $targetdir/$r->{SEED} "; 
  print "$todo\n";

  my $pbs = PBS->new;
  
  $pbs->setScriptName("script_motif_optim_$motif_num.pbs");
  $pbs->setWallTime($walltime);
  
  # now just add commands (that you would normally exexute manually)
  $pbs->addCmd("cd $pwd");
  $pbs->addCmd("echo \"Optimizing $r->{MOTIF}\"");
  $pbs->addCmd($todo);
  $pbs->addCmd("echo \"Done.\"");
  
  if ($submit == 1) {
    # submit the script
    $pbs->submit; 
  } else {
    # alternative: just execute the script, without submitting (useful for debugging)
    $pbs->execute;
  }

  $motif_num ++;
}
