#!/usr/bin/perl

BEGIN{ 

  if ((!$ENV{FIREDIR}) || ($ENV{FIREDIR} eq '')) {
    print "The FIREDIR environment variable is not set. Please set it using export or setenv (see FIRE tutorial online).\n";
    exit;
  }

}

my $cmdline = "fire.pl";
foreach my $r (@ARGV) {
  $cmdline .= " $r";
}

use lib "$ENV{FIREDIR}/SCRIPTS";

my $pwd     = `pwd`; $pwd =~ s/\n//;

my $firedir = $ENV{'FIREDIR'};


my $rootdir        = ".";
my $scriptdir      = "$firedir/SCRIPTS";
my $progdir        = "$firedir/PROGRAMS";       


use Sets;
use PBS;
use Table;
use Getopt::Long;
use strict;


umask 0000;

my $random                 = undef;
my $motiffile_dna          = undef;
my $motiffile_rna          = undef;
my $expfiles               = undef;
my $fa                     = undef;
my $fa_dist                = undef;
my $target_dir             = undef;
my $divbins                = undef;
my $submit                 = 0;
my $quantized              = 1;
my $mbins_dist             = 8;
my $mbins_interval         = 150;
my $add                    = 1;
my $add3                   = undef;
my $add5                   = undef;
my $rna                    = 0;
my $platform               = undef;
my $distonly               = 0;
my $shuffle                = 10000;
my $shuffle_mifind         = -1;
my $shuffle_midist         = 10000;

my $optimslow              = 0;
my $consfile               = undef;
my $removecols_draw        = 0;
my $kmerfile               = undef;
my $rightlabel             = '';
my $leftlabel              = '';
my $seqlen                 = undef;
my $ebins                  = undef;
my $gap                    = undef;
my $maxgocatsize           = 500;
my $walltime               = "20:00:00";
my $debug                  = 0;
my $queue                  = undef;

my $dodef                  = 1;
if (grep(/\-\-dodef=0/, @ARGV) > 0) {
  $dodef = 0;
} 
my $doremovedups           = $dodef;
my $domifind               = $dodef;
my $domioptimize           = $dodef;
my $domisignif             = $dodef;
my $domidist               = $dodef;
my $domisummarize          = $dodef;
my $domicombine            = $dodef;
my $domidrawmatrix         = $dodef;
my $domidrawdist           = $dodef;
my $dogetprofiles          = $dodef;
my $docons                 = $dodef;
my $domotifnames           = $dodef;
my $dogoclusters           = $dodef;
my $domisearch             = undef;
my $dogomotifs             = $dodef;
my $domotifreport          = $dodef;
my $domidrawinteractions   = $dodef;
my $doskipdiscovery        = 0;
my $doalignace             = 0;
my $dobinary               = $dodef;
my $doreduce               = 0;
my $doagglo                = 0;
my $dodrawmotifmaps        = $dodef;
my $doskipoptimization     = 0;
my $dotargets              = $dodef;
my $dohistograms           = $dodef;

my $goindex                = undef;
my $gonames                = undef;

my $nodups                 = 0;

my $donbcopies             = 1;
my $docreateseqs           = 1;
my $fast                   = 1;
my $sortrowsbyphase        = 0;
my $sortmotifsbyphase      = 0;
my $minr                   = 5;
my $minp                   = undef;
my $minz                   = undef;

my $species                = undef;
my $exptype                = "discrete";
my $fastafile_dna          = undef;
my $fastafile_rna          = undef;
my $consfile_dna           = undef;
my $consfile_rna           = undef;
my $shuffle_cons           = 0;
my $rightlabel_dna         = undef;
my $leftlabel_dna          = undef;
my $seqlen_dna             = undef;
my $maxselected            = undef;
my $rightlabel_rna         = undef;
my $leftlabel_rna          = undef;
my $seqlen_rna             = undef;

my $maxfreq                = 0.33;
my $jn                     = 10;
my $jn_t                   = 6;
my $jn_f                   = 3;
my $poscor                 = 1;
my $doallstats             = 1;
my $outlevel               = 2;
my $dorna                  = 1;
my $dodna                  = 1;
my $dodnarna               = 1;

my $fc_fasta1_dna          = undef;
my $fc_fasta2_dna          = undef;
my $fc_nbgenes_dna         = undef;

my $fc_fasta1_rna          = undef;
my $fc_fasta2_rna          = undef;
my $fc_nbgenes_rna         = undef;

my $kmerfile_dna           = undef;
my $kmerfile_rna           = undef;

my $acelibrary_dna         = undef;
my $acelibrary_rna         = undef;

my $micrornas              = undef;

my $jobid_dna              = undef;
my $jobid_rna              = undef;
my $jobid_dnarna           = undef;

my $seqlen_u_dna           = undef;
my $seqlen_d_dna           = undef;

my $oribiasonly            = 1;

my $k                      = undef;
my $suffix                 = undef;
my $checkonly              = 0;
my $dodnacomp              = 1;
my $dornacomp              = 1;
my $mapping_orth           = undef;
my $lp_t_draw              = undef;
my $maxdegeneracy          = undef;


if (@ARGV == 0) {
  die "Usage: [perl] fire.pl --expfiles=FILES --exptype=TXT --species=TXT\n";
}


GetOptions ('expfiles=s'             => \$expfiles,
	    'exptype=s'              => \$exptype,
	    'species=s'              => \$species,	
	    'checkonly=s'            => \$checkonly,
	    'random=s'               => \$random,
	    'motiffile_dna=s'        => \$motiffile_dna,
	    'motiffile_rna=s'        => \$motiffile_rna,	    
	    'consfile_dna=s'         => \$consfile_dna,
	    'consfile_rna=s'         => \$consfile_rna,
	    'kmerfile_dna=s'         => \$kmerfile_dna,
	    'kmerfile_rna=s'         => \$kmerfile_rna,
	    'kmerfile=s'             => \$kmerfile,
	    'fastafile_dna=s'        => \$fastafile_dna,
	    'fastafile_rna=s'        => \$fastafile_rna,
	    'maxselected=s'          => \$maxselected,
	    'k=s'                    => \$k,
	    'gap=s'                  => \$gap,
	    'fa_dist=s'              => \$fa_dist,
	    'walltime=s'             => \$walltime,
	    'target_dir=s'           => \$target_dir,
	    'submit=s'               => \$submit,
	    'rna=s'                  => \$rna,
	    'platform=s'             => \$platform,
	    'distonly=s'             => \$distonly,
	    'shuffle=s'              => \$shuffle,
	    'shuffle_mifind=s'       => \$shuffle_mifind,
	    'mbins_interval=s'       => \$mbins_interval,
	    'shuffle_midist=s'       => \$shuffle_midist,
	    'shuffle_cons=s'         => \$shuffle_cons,
	    'quantized=s'            => \$quantized,
	    'add=s'                  => \$add,
	    'add3=s'                 => \$add3,
	    'add5=s'                 => \$add5,
	    'dodef=s'                => \$dodef,
	    'rootdir=s'              => \$rootdir,
	    'optimslow=s'            => \$optimslow,
	    'removecols_draw=s'      => \$removecols_draw,
	    'nodups=s'               => \$nodups,
	    'seqlen=s'               => \$seqlen,
	    'seqlen_dna=s'           => \$seqlen_dna,
	    'seqlen_rna=s'           => \$seqlen_rna,

	    'leftlabel=s'            => \$leftlabel,
	    'leftlabel_dna=s'        => \$leftlabel_dna,
	    'rightlabel_dna=s'       => \$rightlabel_dna,	    

	    'rightlabel=s'           => \$rightlabel,
	    'leftlabel_rna=s'        => \$leftlabel_rna,
	    'rightlabel_rna=s'       => \$rightlabel_rna,	    

	    'domisearch=s'           => \$domisearch,	    
	    'doremovedups=s'         => \$doremovedups,
	    'domifind=s'             => \$domifind,
	    'domioptimize=s'         => \$domioptimize,
	    'domisignif=s'           => \$domisignif,
	    'domidist=s'             => \$domidist,
	    'domisummarize=s'        => \$domisummarize,
	    'domicombine=s'          => \$domicombine,
	    'domidrawmatrix=s'       => \$domidrawmatrix,
	    'domidrawinteractions=s' => \$domidrawinteractions,
	    'dogoclusters=s'         => \$dogoclusters,
	    'dogomotifs=s'           => \$dogomotifs,
	    'dogetprofiles=s'        => \$dogetprofiles,
	    'doagglo=s'              => \$doagglo,
	    'dodna=s'                => \$dodna,
	    'dorna=s'                => \$dorna,
	    'dodnarna=s'             => \$dodnarna,
	    'docons=s'               => \$docons,
	    'domotifnames=s'         => \$domotifnames,	    
	    'doskipdiscovery=s'      => \$doskipdiscovery,
	    'domidrawdist=s'         => \$domidrawdist,
	    'docreateseqs=s'         => \$docreateseqs,
	    'doalignace=s'           => \$doalignace,
	    'doreduce=s'             => \$doreduce,
	    'dodrawmotifmaps=s'      => \$dodrawmotifmaps,
	    'doskipoptimization=s'   => \$doskipoptimization,
	    'dotargets=s'            => \$dotargets,
	    'dohistograms=s'         => \$dohistograms,
	    'sortrowsbyphase=s'      => \$sortrowsbyphase,
	    'fast=s'                 => \$fast,
	    'minr=s'                 => \$minr,
	    'minp=s'                 => \$minp,
	    'minz=s'                 => \$minz,
	    'divbins=s'              => \$divbins,
	    'maxfreq=s'              => \$maxfreq,
	    'domotifreport=s'        => \$domotifreport,
	    'jn=s'                   => \$jn,
	    'jn_t=s'                 => \$jn_t,
	    'jn_f=s'                 => \$jn_f,
	    'poscor=s'               => \$poscor,
	    'doallstats=s'           => \$doallstats,
	    'sortmotifsbyphase=s'    => \$sortmotifsbyphase,
	    'mbins_dist=s'           => \$mbins_dist,
	    'oribiasonly=s'          => \$oribiasonly,
	    'ebins=s'                => \$ebins,
	    'suffix=s'               => \$suffix,
	    'dodnacomp=s'            => \$dodnacomp,
	    'dornacomp=s'            => \$dornacomp,
	    'debug=s'                => \$debug,
	    'queue=s'                => \$queue,
	    'maxgocatsize=s'         => \$maxgocatsize,
	    'seqlen_u_dna=s'         => \$seqlen_u_dna,
	    'seqlen_d_dna=s'         => \$seqlen_d_dna,
	    'maxdegeneracy=s'        => \$maxdegeneracy,
	    'lp_t_draw=s'            => \$lp_t_draw);

if (!defined($expfiles)) { # && !defined($expfile)) {
  die("Please input an expression file (--expfiles=FILE).\n");
}

if (defined($domisearch) && ($domisearch == 0)) {

  $doremovedups   = 0;
  $domifind       = 0;
  $domioptimize   = 0;
  $domisignif     = 0;
  $domidist       = 0;

} elsif (defined($domisearch) && ($domisearch == 1))  {
  
  $doremovedups   = 1;
  $domifind       = 1;
  $domioptimize   = 1;
  $domisignif     = 1;
  $domidist       = 1;

}

if ($doskipdiscovery == 1) {
  $domifind       = 0;
  $domioptimize   = 0;
}

if ($doskipoptimization == 1) {
  $domioptimize   = 0;
}

if (defined($exptype)) {
  if ($exptype eq "continuous") {
    $quantized    = 0;
  } elsif ($exptype eq "discrete") {
    $quantized    = 1;
  } else {
    die "--exptype can be either discrete or continuous\n";
  }
}

if (($quantized == 0) && (!defined($removecols_draw))) {
  $removecols_draw = 0;
}




if (defined($species)) {

  #
  #  read species file
  #
  my $species_data = readSpeciesData($species);
  
  if (-e $species_data->{"mapping_orth"}) {
    $mapping_orth = $species_data->{"mapping_orth"};
  }
  
  if (!defined($fastafile_dna)) {
    $fastafile_dna = $species_data->{"fastafile_dna"};
  }  
  
  if (!defined($fastafile_rna)) {
    $fastafile_rna = $species_data->{"fastafile_rna"};
  }

  if (!defined($kmerfile_dna)) {
    $kmerfile_dna = $species_data->{"kmerfile_dna"};
  }  
  
  if (!defined($kmerfile_rna)) {
    $kmerfile_rna = $species_data->{"kmerfile_rna"};
  }

  if ($rightlabel_dna eq "") {
    $rightlabel_dna = $species_data->{"rightlabel_dna"};
  }
  
  if ($leftlabel_dna eq "") {
    $leftlabel_dna  = $species_data->{"leftlabel_dna"};
  }

  if (!defined($seqlen_dna)) {
    $seqlen_dna = $species_data->{"seqlen_dna"};
  }


  if (defined($seqlen_u_dna) && defined($seqlen_d_dna)) {
    $fastafile_dna    =  $species_data->{"fastafile_dna_p"} . "$seqlen_u_dna\_$seqlen_d_dna.fa";
    die "$fastafile_dna does not exist, sorry.\n" if (! -e $fastafile_dna);
    $seqlen_dna       =  $species_data->{"seqlen_dna"};
    $leftlabel_dna    = -$species_data->{"seqlen_u_dna"};
    $rightlabel_dna   =  $species_data->{"seqlen_d_dna"};
  }

  if ($rightlabel_rna eq "") {
    $rightlabel_rna = $species_data->{"rightlabel_rna"};
  }
  
  if ($leftlabel_rna eq "") {
    $leftlabel_rna  = $species_data->{"leftlabel_rna"};
  }
  
  if (!defined($seqlen_rna)) {
    $seqlen_rna     = $species_data->{"seqlen_rna"};
  }
  
  $goindex = $species_data->{"goindex"} if  (!defined($goindex));
  $gonames = $species_data->{"gonames"} if  (!defined($gonames));

  $fc_fasta1_dna  = $species_data->{"fc_fasta1_dna"} if (defined($species_data->{"fc_fasta1_dna"}));
  $fc_fasta2_dna  = $species_data->{"fc_fasta2_dna"} if (defined($species_data->{"fc_fasta2_dna"}));
  $fc_nbgenes_dna = $species_data->{"fc_nbgenes_dna"} if (defined($species_data->{"fc_nbgenes_dna"}));

  $fc_fasta1_rna  = $species_data->{"fc_fasta1_rna"} if (defined($species_data->{"fc_fasta1_rna"}));
  $fc_fasta2_rna  = $species_data->{"fc_fasta2_rna"} if (defined($species_data->{"fc_fasta2_rna"}));
  $fc_nbgenes_rna = $species_data->{"fc_nbgenes_rna"} if (defined($species_data->{"fc_nbgenes_rna"}));

  if (defined($species_data->{"acelibrary_dna"})) {
    $acelibrary_dna = $species_data->{"acelibrary_dna"};
  }
  
  if (defined($species_data->{"acelibrary_rna"})) {
    $acelibrary_rna = $species_data->{"acelibrary_rna"};
  }

  if (defined($species_data->{"micrornas"})) {
    $micrornas = $species_data->{"micrornas"};
  }
  
  
}


my $time = Sets::getNiceDateTime(1);

if (defined($suffix)) {

  my $a_ref_files = Sets::getFiles($expfiles);
  if (@$a_ref_files > 1) {
    die "Cannot run -suffix on multiple files.\n";
  }
  
  system("cp $expfiles $expfiles.$suffix");

  $expfiles .= ".$suffix";

}

if (defined($target_dir)) {

  # create target_dir if needed
  if (!-e $target_dir) {
    mkdir $target_dir;
  }
  
  # copy expfiles
  system("cp $expfiles $target_dir/$expfiles");

  # change expfile
  $expfiles = "$target_dir/$expfiles";

}





if (defined($random)) {
  
  my $a_ref_files = Sets::getFiles($expfiles);
  if (@$a_ref_files > 1) {
    die "Cannot run -random on multiple files.\n";
  }
  
  my $target_dir = "$expfiles\_RANDOM"; 
  system("rm -Rf $target_dir");
  
  mkdir $target_dir if (! -e $target_dir);
  my $expfile_file = Sets::filename($expfiles);

  for (my $i=0; $i<$random; $i++) {
    my $cmd = "perl $scriptdir/shuffle_column.pl $expfiles > $target_dir/$expfile_file.$i.txt";
    print "$cmd\n";
    system($cmd);
  }

  $expfiles .= "\_RANDOM/*.txt";

  $domicombine            = 0;
  $domidrawmatrix         = 0;
  $domidrawdist           = 0;
  $dogetprofiles          = 0;
  $docons                 = 0;
  $domotifnames           = 0;
  $dogoclusters           = 0;
  $dogomotifs             = 0;
  $domotifreport          = 0;
  $domidrawinteractions   = 0;

}


#
#  go over the expression files
# 
my $a_ref_files = Sets::getFiles($expfiles);
print "FF=$expfiles\n";

foreach my $expfile (@$a_ref_files) {

  my $todo = 'perl -pi -e "s/[\r\ ]//g" ' . $expfile;
  system($todo);

  if (&check_input_file($expfile, $exptype) == 0) {
    die "Please correct input expression file.\n";
  }
  
  if ($checkonly == 1) {
    next;
  }

  my $expfile_file = Sets::filename($expfile);
  
  $target_dir = "$expfile\_FIRE";

  if (! -e $target_dir) {
    mkdir $target_dir;
  }

  #
  # save command line to _FIRE
  # 
  open OUTC, ">$target_dir/cmdline.txt" or print "Cannot open $target_dir/cmdline.txt";
  print OUTC "$cmdline\n";
  close OUTC;

  #
  #
  # **********         DNA ANALYSIS        ***********
  #
  #

  my $target_dir_dna = "$target_dir/DNA";
  if (! -e $target_dir_dna) {
    mkdir $target_dir_dna; 
  }

  my $outdir_seq_dna          = "$target_dir_dna/sequences";
  my $expfile_nodups_dna      = "$target_dir_dna/$expfile_file";
  my $expfile_q_dna           = "$expfile_nodups_dna.quantized";
  my $quantized_expfile_dna   = $expfile_nodups_dna;
  if ($quantized == 0) {
    $quantized_expfile_dna = $expfile_q_dna;
  }

  my $expfile_b_dna           = "$expfile_nodups_dna.binary";
  my $expfile_b_summary_dna   = "$expfile_nodups_dna.summary.binary";
  my $seedfile_dna            = "$expfile_nodups_dna.seeds";
  my $optimfile_dna           = "$expfile_nodups_dna.optim";
  my $distfile_dna            = "$expfile_nodups_dna.dist";
  my $signiffile_dna          = "$expfile_nodups_dna.signif";
  my $motifrepfile_dna        = "$expfile_nodups_dna.signif.motifs.rep";
  my $summaryfile_dna         = "$expfile_nodups_dna.summary";
  my $densityfile_dna         = "$expfile_nodups_dna.densities";
  my $epsdensityfile_dna      = "$expfile_nodups_dna.densities.eps";
  my $matrixfile_dna          = "$expfile_nodups_dna.matrix";
  my $fullmatrixfile_dna      = "$expfile_nodups_dna.fullmatrix";  
  my $columnsfile_dna         = "$expfile_nodups_dna.columns";
  my $signifcolumnsfile_dna   = "$expfile_nodups_dna.signifcolumns";
  my $namesfile_dna           = "$expfile_nodups_dna.motifnames";
  my $clusterfile_dna         = "$expfile_nodups_dna.clusters";
  my $mimatrixfile_dna        = "$expfile_nodups_dna.mimatrix";
  my $fullmimatrixfile_dna    = "$expfile_nodups_dna.fullmimatrix";
  my $profiles_dna            = "$expfile_nodups_dna.profiles";
  my $gofile_dna              = "$expfile_nodups_dna.GO";
  my $gofile_full_dna         = "$expfile_nodups_dna.GO.full";
  my $gofile_pairs_dna        = "$expfile_nodups_dna.GOpairs";
  my $gofile_pairs_full_dna   = "$expfile_nodups_dna.GOpairs.full";
  my $consfile_dna            = "$expfile_nodups_dna.cons";
  my $gomofile_dna            = "$expfile_nodups_dna.GOmotifs";
  my $gomofile_full_dna       = "$expfile_nodups_dna.GOmotifs.full";
  my $outdistmatrix_dna       = "$expfile_nodups_dna.distmatrix";
  my $outoriematrix_dna       = "$expfile_nodups_dna.oriematrix";
  my $orth_profiles_dna       = "$expfile_nodups_dna.profiles_orth";
  my $motifreport_dna         = "$expfile_nodups_dna.motifreport";
  my $targets_dna             = "$expfile_nodups_dna.targets";
  my $histogramdir_dna        = "$expfile_nodups_dna.histograms_OUT";
 
 if ($dodna == 1) {
  
    #
    # make sure sequence file exists
    #
    die "No DNA sequence data ($fastafile_dna does not exist). Perhaps you need to download the relevant 
species-specific data file from http://tavazoielab.princeton.edu/FIRE.\n" if (! -e $fastafile_dna);
    
    #
    # start script
    #
    my $pbs = PBS->new;
    
    $pbs->setPlatform($platform) if (defined($platform));
    $pbs->setQueue($queue)       if (defined($queue));

    #$pbs->setMemory("4096Mb");

    $pbs->setWallTime($walltime);

    if ($platform eq 'tcluster') {
      $pbs->addCmd("setenv FIREDIR $firedir");	      
      $pbs->addCmd("setenv LD_LIBRARY_PATH$ENV{FIREDIR}/modules/lib");      
    } else {
      $pbs->addCmd("export FIREDIR=$firedir");	
      $pbs->addCmd("export LD_LIBRARY_PATH=$ENV{FIREDIR}/modules/lib");
    }

    $pbs->addCmd("cd $pwd");

    

    #if (defined($platform) && ($platform eq "fafner")) {
    #  $pbs->addCmd("export DYLD_LIBRARY_PATH=/Genomics/fafner/grid/users/elemento/usr/lib");
    #}
    
 

    $pbs->setScriptName("$expfile_nodups_dna.script");

    $pbs->addCmd("date");
    $pbs->addCmd("echo \"DNA, remove duplicates, create $expfile_nodups_dna\"");

    if ($nodups == 1) {
      system("cp $expfile $expfile_nodups_dna");

    } elsif  ($doremovedups == 1) {

      my %PARAMS = ("expfile"       => $expfile,
		    "quantized"     => $quantized, 
		    "fastafile"     => $fastafile_dna,
		    "dupfile"       => undef,
		    "ebins"         => $ebins,	
		    "divbins"       => $divbins,	    
		    "outfile"       => $expfile_nodups_dna);

      my $cmd = &get_cmd_removedups(\%PARAMS); 
      $pbs->addCmd($cmd);
    }

    #
    # quantize if necessary
    #
    if ($quantized == 0) {
      $pbs->addCmd("echo \"DNA, quantizing.\"");    
      my $cmd = "perl $scriptdir/quantize_expression_vector.pl -expfile $expfile_nodups_dna -outfile $expfile_q_dna ";	
      if (defined($ebins)) {
	$cmd .= " -ebins $ebins ";
      }
      if (defined($divbins)) {
	$cmd .= " -divbins $divbins ";
      }
      $pbs->addCmd($cmd);
    }
    
    
    #
    # seeds
    #
    if ($domifind == 1) { 
      $pbs->addCmd("echo \"DNA, Step 1: seed discovery.\""); 

      my %PARAMS = ("expfile"     => $expfile_nodups_dna,
		    "quantized"   => $quantized,
		    "fastafile"   => $fastafile_dna, 
		    "shuffle"     => $shuffle_mifind, 
		    "kmerfile"    => $kmerfile,
		    "fast"        => $fast,
		    "rna"         => 0,
		    "maxselected" => $maxselected,
		    "gap"         => $gap,
		    "k"           => $k,
		    "ebins"       => $ebins,
		    "divbins"     => $divbins,
		    "seedfile"    => $seedfile_dna);
      my $cmd = &get_cmd_mifind(\%PARAMS);
      $pbs->addCmd($cmd);
    }
  
    # if no seeds, exit here
    if ($doskipdiscovery == 0) {
      $pbs->addCmd("CNTM=`perl $scriptdir/count_lines.pl $seedfile_dna`; if [ \$CNTM = 0 ]; then echo \"No DNA motifs. Exiting now.\"; exit; fi");      
    }

    if ($domioptimize == 1) {

      $pbs->addCmd("echo \"DNA, Step 2: seed optimization.\""); 
    
      my %PARAMS = ("expfile"     => $expfile_nodups_dna, 
		    "quantized"   => $quantized, 
		    "fastafile"   => $fastafile_dna, 
		    "rna"         => 0, 
		    "gap"         => $gap,
		    "seedfile"    => $seedfile_dna, 
		    "shuffle"     => $shuffle, 
		    "minr"        => $minr, 
		    "add"         => $add,
		    "maxfreq"     => $maxfreq, 
		    "optimslow"   => $optimslow,
		    "maxdegeneracy" => $maxdegeneracy,
		    "ebins"       => $ebins, 
		    "divbins"     => $divbins,
		    "optimfile"   => $optimfile_dna);
      my $cmd = &get_cmd_mioptimize(\%PARAMS);
      $pbs->addCmd($cmd);
    }


    if ($doskipdiscovery == 1) {
      $pbs->addCmd("echo \"DNA, Step 2: skip seed optimization, use -motiffile_dna option instead.\"");
      die "Please provide -motiffile_dna\n" if (! -e $motiffile_dna);
      $pbs->addCmd("perl $scriptdir/copy_seeds_to_optim.pl $motiffile_dna $optimfile_dna");
    }


    if ($domisignif == 1) {
      $pbs->addCmd("echo \"DNA, Step 3: evaluation of motif significances.\""); 

      my %PARAMS = ("expfile"     => $expfile_nodups_dna, 
		    "quantized"   => $quantized, 
		    "fastafile"   => $fastafile_dna, 
		    "rna"         => 0,
		    "optimfile"   => $optimfile_dna, 
		    "shuffle"     => $shuffle, 
		    "jn"          => $jn, 
		    "jn_t"        => 0, 
		    "jn_f",       => $jn_f,
		    "ebins"       => $ebins,
		    "divbins"     => $divbins,
		    "signiffile"  => $signiffile_dna);
      my $cmd = &get_cmd_misignif(\%PARAMS);
      $pbs->addCmd($cmd);
    
    }
    
    if ($dobinary == 1) {

      my $cmd = "perl $scriptdir/mi_create_overep_expfile.pl -repfile $motifrepfile_dna -expfile $quantized_expfile_dna -quantized $quantized -outexpfile $expfile_b_dna";
      $pbs->addCmd("echo \"DNA, Step 3.5: creating binary expression profiles.\"");
      $pbs->addCmd($cmd);

    }


    if ($domidist == 1) {
      $pbs->addCmd("echo \"DNA, Step 4: discovery of distance constraints.\""); 
    
      my %PARAMS = ("expfile_b"      => $expfile_b_dna, 
		    "expfile"        => $expfile_nodups_dna, 
		    "quantized"      => $quantized, 
		    "fastafile"      => $fastafile_dna, 
		    "rna"            => 0, 
		    "optimfile"      => $optimfile_dna, 
		    "shuffle"        => $shuffle, 
		    "mbins_interval" => $mbins_interval, 
		    "donbcopies"     => $donbcopies, 
		    "ebins"          => $ebins,
		    "divbins"        => $divbins,
		    "distfile"       => $distfile_dna);
      my $cmd = &get_cmd_midist(\%PARAMS);
      $pbs->addCmd($cmd);
    }


    if (($dogoclusters == 1) && (defined($goindex))) {
      
      $pbs->addCmd("echo \"DNA, Step 6.5: GO analysis.\"");    

      my $cmd = "perl $scriptdir/clusters_go_enrichment.pl --clusters=$quantized_expfile_dna --goindex=$goindex --gonames=$gonames --N=-1 --outfile=$gofile_full_dna --outgo=$gofile_dna --maxcatsize=$maxgocatsize ";

      #
      #  fafner code to fix malfunctioning Hypergeom::cumhyper on fafner only
      #
      if (defined($platform) && ($platform eq "fafner")) {
	$cmd .= " --usemodule=0 ";
      }

      $pbs->addCmd($cmd);    

    }


    if (($docons == 1) && (defined($fc_fasta1_dna))) {

      $pbs->addCmd("echo \"DNA, Step 5.5: conservation.\"");    

      my $cmd = "perl $scriptdir/evaluate_motifs_conservation.pl -summaryfile $optimfile_dna -fasta1 $fc_fasta1_dna -fasta2 $fc_fasta2_dna -nbgenes $fc_nbgenes_dna -kmerfile $kmerfile_dna -outfile $consfile_dna ";
      if ($shuffle_cons == 1) {
	$cmd .= " -shuffle 1";
      }
      $pbs->addCmd($cmd);
    }
    
  
    if ($domisummarize == 1) {
      $pbs->addCmd("echo \"DNA, Step 5: summarize information.\"");

      my %PARAMS = ("expfile_dna"        => $expfile_nodups_dna,
		    "fastafile_dna"      => $fastafile_dna,
		    "optimfile_dna"      => $optimfile_dna,
		    "signiffile_dna"     => $signiffile_dna,
		    "distfile_dna"       => $distfile_dna,
		    "consfile_dna"       => (defined($fc_fasta1_dna)?$consfile_dna:undef),
		    "rna"                => 0,
		    "quantized"          => $quantized,
		    "minrobustness"      => $jn_t,
                    "maxrobustness"      => $jn,
		    "removecols_draw"    => $removecols_draw,
		    "sortrowsbyphase"    => $sortmotifsbyphase,
		    "summaryfile"        => $summaryfile_dna,
		    "matrixfile"         => $matrixfile_dna,
		    "fullmatrixfile"     => $fullmatrixfile_dna,
		    "densityfile"        => $densityfile_dna,
		    "columnsfile"        => $columnsfile_dna,
		    "signifcolumnsfile"  => $signifcolumnsfile_dna,
		    "outlevel"           => $outlevel,
		    "minp"               => $minp,
		    "minz"               => $minz,		    
		    "namesfile"          => $namesfile_dna,
		    
		    "expfile_b_dna"      => $expfile_b_dna,
		    "outexpfile_b"       => $expfile_b_summary_dna
		   );

      my $cmd = &get_cmd_misummarize(\%PARAMS);
      $pbs->addCmd($cmd);
  
    }
  
    if ($dogetprofiles == 1) {

      $pbs->addCmd("echo \"DNA, Step 5.6: generate profiles for motifs in summary file.\"");
      my $cmd = "perl $scriptdir/generate_motif_profiles.pl -fastafile $fastafile_dna -summaryfile $summaryfile_dna -outfile $profiles_dna";
      $pbs->addCmd($cmd);
      
      #
      #  generate profiles for orthologs
      #
      if (defined($fc_fasta2_dna)) {
	$pbs->addCmd("echo \"DNA, Step 5.6 bis: Generate profiles for $fc_fasta2_dna.\"");
	my $cmd = "perl $scriptdir/generate_motif_profiles.pl -fastafile $fc_fasta2_dna -summaryfile $summaryfile_dna -outfile $orth_profiles_dna";
	$pbs->addCmd($cmd);
      }
    }
    
    
    if ($dotargets == 1) {
      
      $pbs->addCmd("echo \"DNA, Determining target gene sets.\"");
      my $cmd = "perl $scriptdir/get_motif_targets.pl --matrixfile=$fullmatrixfile_dna --profiles=$profiles_dna --expfile=$quantized_expfile_dna --outfile=$targets_dna";
      $pbs->addCmd($cmd);

    }


    if ($dohistograms == 1) {

      $pbs->addCmd("echo \"DNA, Drawing position histograms.\"");
      my $cmd = " perl $scriptdir/draw_position_histogram.pl  --profiles=$profiles_dna --expfile=$quantized_expfile_dna --outdir=$histogramdir_dna --targetfile=$targets_dna --rna=0 --seqlen=$seqlen_dna --bins=10 ";
      $pbs->addCmd($cmd);

    }


    if (($dodnacomp == 1) && ($domicombine == 1)) {


      $pbs->addCmd("echo \"DNA, Step 6: cluster motifs.\"");
    
      my %PARAMS = ("expfile"               => $expfile_b_summary_dna, 
		    "quantized"             => 1,
		    "fastafile_dna"         => $fastafile_dna, 
		    "summaryfile"           => $summaryfile_dna,
		    "poscor"                => $poscor, 
		    "doallstats"            => $doallstats,
		    "clusterfile"           => $clusterfile_dna,
		    "mimatrixfile"          => $mimatrixfile_dna, 
		    "fullmimatrixfile"      => $fullmimatrixfile_dna);

      my $cmd = &get_cmd_micombine(\%PARAMS);
      $pbs->addCmd($cmd);

    }

    #
    # agglomerative clustering of motifs based on p-values ... 
    #   will overwrite .clusters file (experimental)
    #
    if ($doagglo == 1) {
      
      my $nbclusters = 4;
      my $todo = "perl $scriptdir/mi_cluster_motifs.pl -nbclusters $nbclusters -matrixfile $matrixfile_dna -clusterfile $clusterfile_dna";
      $pbs->addCmd($todo);

    }


    if ($dogomotifs == 1) {
       
      $pbs->addCmd("echo \"DNA, Step 5.6: generate genome-wide GO enrichments for motifs in summary file.\"");

      my $cmd = "perl $scriptdir/motif_single_go_enrichment.pl --summaryfile=$summaryfile_dna --profiles=$profiles_dna --goindex=$goindex --gonames=$gonames --outgo=$gomofile_dna --outfile=$gomofile_full_dna --expfile=$quantized_expfile_dna --maxcatsize=$maxgocatsize ";
      
      #
      #  fafner code to fix malfunctioning Hypergeom::cumhyper on fafner only
      #
      if (defined($platform) && ($platform eq "fafner")) {
	$cmd .= " --usemodule=0 ";
      }


      if (($dodnacomp == 1) && ($domicombine == 1)) {
	$cmd .= " --clusterfile=$clusterfile_dna";
      }
      
      $pbs->addCmd($cmd);

      
    }
  
    if ($dodrawmotifmaps == 1) {
      #
      # draw interaction motif maps
      #
      $pbs->addCmd("echo \"DNA, draw motif maps for co-localizing pairs\"");
      mkdir "$mimatrixfile_dna\_OUT" if (! -e "$mimatrixfile_dna\_OUT");

      my $cmd = "perl $scriptdir/mi_draw_all_pair_motif_maps.pl --expfile=$quantized_expfile_dna --profiles=$profiles_dna --seqlen=$seqlen_dna --rna=0 --summaryfile=$summaryfile_dna --leftlabel=$leftlabel_dna --rightlabel=$rightlabel_dna --rootdir=$rootdir --fullmatrixfile=$fullmatrixfile_dna ";
      
      if (($domicombine >= 1) && ($dodnacomp == 1)) {
	$cmd .= " --mimatrixfile=$mimatrixfile_dna "; 
      }

      $pbs->addCmd($cmd);
    }

  

    
    if ($domotifnames == 1) {
      
      if (defined($acelibrary_dna)) {
	my $namesfile_dna_tmp = "$namesfile_dna.tmp";
	$pbs->addCmd("echo \"DNA, Step 6.7: identifying motifs.\"");
	my $cmd = "perl $scriptdir/annotate_motifs_using_ace_library.pl --acefiles=\"$acelibrary_dna\" --summaryfile=$summaryfile_dna --namefile=$namesfile_dna --rna=0 --outfile=$namesfile_dna_tmp";	
	$pbs->addCmd($cmd);
	$pbs->addCmd("cp $namesfile_dna $namesfile_dna.old");
	$pbs->addCmd("cp $namesfile_dna.tmp $namesfile_dna");
      }
      
    }

    if (($dodnacomp == 1) && ($domidrawmatrix == 1)) {
      $pbs->addCmd("echo \"DNA, Step 7: draw matrix figure.\"");

      my %PARAMS = ("matrixfile"      => $matrixfile_dna,
		    "summaryfile"     => $summaryfile_dna,
		    "columnsfile"     => $columnsfile_dna,
		    "expfile"         => $quantized_expfile_dna,
		    "gofile"          => (defined($goindex)?$gofile_dna:undef),
		    "ps2pdf"          => 1,
		    "every"           => 1,
		    "lp_t_draw"       => $lp_t_draw,
		    "quantized"       => $quantized,
		    "colmap"          => "$scriptdir/HEATMAPS/cmap2.txt",
		    "namesfile"       => $namesfile_dna);
      
      $PARAMS{"clusterfile"} = $clusterfile_dna if (($sortmotifsbyphase == 0) && ($domicombine >= 1)) ;


      my $cmd  = &get_cmd_midrawmatrix(\%PARAMS);
      $pbs->addCmd($cmd);

      #
      # draw density figure
      #
      %PARAMS    = ("matrixfile"      => $matrixfile_dna,
                    "densityfile"     => $densityfile_dna,
		    "summaryfile"     => $summaryfile_dna,
		    "columnsfile"     => $columnsfile_dna,
		    "expfile"         => $quantized_expfile_dna,
		    "gofile"          => (defined($goindex)?$gofile_dna:undef),		    
		    "ps2pdf"          => 1,
		    "every"           => 1,
		    "quantized"       => $quantized,
		    "redoweblogo"     => 0,
		    "colmap"          => "$scriptdir/HEATMAPS/cmap_dens.txt",
		    "outeps"          => $epsdensityfile_dna,
		    "namesfile"       => $namesfile_dna);
      $PARAMS{"clusterfile"} = $clusterfile_dna if (($sortmotifsbyphase == 0) && ($domicombine >= 1)) ;

      $cmd  = &get_cmd_midrawmatrix(\%PARAMS);
      $pbs->addCmd($cmd);

    }

    
    if ($domidrawdist == 1) {
      
      $pbs->addCmd("echo \"DNA, Step 7b: some distance stuff.\"");
      
      # draw position/orientation biases maps when bias found
      
      my $todo = "perl  $scriptdir/mi_create_dist_pvalues.pl -expfile $quantized_expfile_dna -profiles $profiles_dna -distrepfile $distfile_dna.rep -summaryfile $summaryfile_dna -motifdir $summaryfile_dna\_OUT -quantized $quantized -seqlen $seqlen_dna -rightlabel $rightlabel_dna -leftlabel $leftlabel_dna -outdistmatrix $outdistmatrix_dna -outoriematrix $outoriematrix_dna ";
      $pbs->addCmd($todo);
      
    }
    

    if (($dodnacomp == 1) && ($domidrawinteractions == 1)) {

      my %PARAMS= ("summaryfile"       => $summaryfile_dna,
		   "fullmimatrixfile"  => $fullmimatrixfile_dna,
		   "mimatrixfile"      => $mimatrixfile_dna,
		   "motifdir"          => "$summaryfile_dna\_OUT");

      $PARAMS{"clusterfile"} = $clusterfile_dna if (($sortmotifsbyphase == 0) && ($domicombine >= 1));      
      if ($sortmotifsbyphase == 1) {
	$PARAMS{"orderfile"} = $matrixfile_dna;
      }
      my $cmd = &get_cmd_midrawinteractions(\%PARAMS);
      $pbs->addCmd($cmd);
    }


    if ($domotifreport == 1) {

      $pbs->addCmd("echo \"DNA, Step 10: motif report.\"");

      my $cmd = "perl $scriptdir/generate_motif_occurrence_summary.pl -expfile $expfile_nodups_dna -summaryfile $summaryfile_dna -profiles $profiles_dna -rna 0  -pvmatrix $matrixfile_dna -columnfile $columnsfile_dna -distmatrix $outdistmatrix_dna -oriematrix $outoriematrix_dna -quantized $quantized ";
      $cmd .= " -gomotifs $gomofile_dna -goindex $goindex " if (defined($goindex));
      $cmd .= " -profiles_orth $orth_profiles_dna " if (defined($fc_fasta2_dna));
      $cmd .= " -mapping_orth $mapping_orth " if (defined($mapping_orth));

      if ($quantized == 0) {
	$cmd .= " -expfile_q $expfile_q_dna ";
      }

      $cmd .= " > $motifreport_dna";
      $pbs->addCmd($cmd);
  
    }

    #
    # OBSOLETE : used for publication only
    #

    #if ($doalignace == 1) {      
    #  my $cmd = "perl $scriptdir/clusters_run_alignace.pl --clusters=$expfile_nodups_dna --fastafile=$fastafile_dna > $expfile_nodups_dna.AlignACE";
    #  $pbs->addCmd($cmd);
    #}

    #    if ($doreduce == 1) {    
    #      my $td = "$expfile_nodups_dna.REDUCE_OUT";
    #      mkdir $td if (! -e $td);
    #      my $cmd = "MATRIXREDUCE/Linux/MatrixREDUCE/bin/MatrixREDUCE -sequence=$fastafile_dna -expression=$expfile_nodups_dna  -output=$td -runlog=stdout -p_value=0.0001 -dyad_length=3 -min_gap=1 -max_gap=1 -flank=1 > $expfile_nodups_dna.REDUCE";  #-single_strand
    #      $pbs->addCmd($cmd);
    #    }
    
    #
    # END OBSOLETE STUFF
    #

    
    $pbs->addCmd("date");

    if ($debug == 1) {
       $pbs->print;
    }

    if ($submit == 0) {
      $pbs->execute;
    } else {
      $jobid_dna = $pbs->submit; print "Submitted job $jobid_dna.\n";
    }
    
  }
  #
  #
  #    ***** end DNA *****
  #
  #
  
  
  
  #
  #
  #    *****      RNA         **** 
  #
  #
  
  

  #
  # remove RNA duplicates
  #
  
  my $target_dir_rna = "$target_dir/RNA";
  if (! -e $target_dir_rna) {
    mkdir $target_dir_rna; 
  }

  my $outdir_seq_rna          = "$target_dir_rna/sequences";
  my $expfile_nodups_rna      = "$target_dir_rna/$expfile_file";
  my $expfile_q_rna           = "$expfile_nodups_rna.quantized";
  my $quantized_expfile_rna   = $expfile_nodups_rna;
  if ($quantized == 0) {
    $quantized_expfile_rna = $expfile_q_rna;
  }

  my $expfile_b_rna           = "$expfile_nodups_rna.binary";
  my $expfile_b_summary_rna   = "$expfile_nodups_rna.summary.binary";
  my $seedfile_rna            = "$expfile_nodups_rna.seeds";
  my $optimfile_rna           = "$expfile_nodups_rna.optim";
  my $distfile_rna            = "$expfile_nodups_rna.dist";
  my $signiffile_rna          = "$expfile_nodups_rna.signif";
  my $motifrepfile_rna        = "$expfile_nodups_rna.signif.motifs.rep";
  my $densityfile_rna         = "$expfile_nodups_rna.densities";
  my $epsdensityfile_rna      = "$expfile_nodups_rna.densities.eps";
  my $summaryfile_rna         = "$expfile_nodups_rna.summary";
  my $matrixfile_rna          = "$expfile_nodups_rna.matrix";
  my $fullmatrixfile_rna      = "$expfile_nodups_rna.fullmatrix";  
  my $columnsfile_rna         = "$expfile_nodups_rna.columns";
  my $signifcolumnsfile_rna   = "$expfile_nodups_rna.signifcolumns";
  my $namesfile_rna           = "$expfile_nodups_rna.motifnames";
  my $clusterfile_rna         = "$expfile_nodups_rna.clusters";
  my $mimatrixfile_rna        = "$expfile_nodups_rna.mimatrix";
  my $fullmimatrixfile_rna    = "$expfile_nodups_rna.fullmimatrix";
  my $profiles_rna            = "$expfile_nodups_rna.profiles";
  my $gofile_rna              = "$expfile_nodups_rna.GO";
  my $gofile_full_rna         = "$expfile_nodups_rna.GO.full";
  my $gofile_pairs_rna        = "$expfile_nodups_rna.GOpairs";
  my $gofile_pairs_full_rna   = "$expfile_nodups_rna.GOpairs.full";
  my $consfile_rna            = "$expfile_nodups_rna.cons";
  my $gomofile_rna            = "$expfile_nodups_rna.GOmotifs";
  my $gomofile_full_rna       = "$expfile_nodups_rna.GOmotifs.full";  
  my $outdistmatrix_rna       = "$expfile_nodups_rna.distmatrix";
  my $outoriematrix_rna       = "$expfile_nodups_rna.oriematrix";
  my $orth_profiles_rna       = "$expfile_nodups_rna.profiles_orth";
  my $motifreport_rna         = "$expfile_nodups_rna.motifreport";
  my $targets_rna             = "$expfile_nodups_rna.targets";
  my $histogramdir_rna        = "$expfile_nodups_rna.histograms_OUT";

  if ($dorna == 1) {
    
  
    #
    # make sure sequence file exists
    #
    die "No RNA sequence data. Perhaps you need to download the relevant 
species-specific data file from http://tavazoielab.princeton.edu/FIRE.\n" if (! -e $fastafile_rna);
  
    #
    # start script
    #
    my $pbs = PBS->new;

    $pbs->setWallTime($walltime);

    $pbs->setQueue($queue)       if (defined($queue));
    $pbs->setPlatform($platform) if (defined($platform));

    $pbs->addCmd("cd $pwd");

    if ($platform eq 'tcluster') {
      $pbs->addCmd("setenv FIREDIR $firedir");	      
      $pbs->addCmd("setenv LD_LIBRARY_PATH$ENV{FIREDIR}/modules/lib");      
    } else {
      $pbs->addCmd("export FIREDIR=$firedir");	
      $pbs->addCmd("export LD_LIBRARY_PATH=$ENV{FIREDIR}/modules/lib");
    }

    #if (defined($platform) && ($platform eq "fafner")) {
    #  $pbs->addCmd("export DYLD_LIBRARY_PATH=/Genomics/fafner/grid/users/elemento/usr/lib");
    #}
    
  
    $pbs->setScriptName("$expfile_nodups_rna.script");

    $pbs->addCmd("date");
    $pbs->addCmd("echo \"RNA, Processing $expfile_nodups_rna\"");


    if ($nodups == 1) {
      system("cp $expfile $expfile_nodups_rna");

    } elsif ($doremovedups == 1) {
      my %PARAMS = ("expfile"       => $expfile,
		    "quantized"     => $quantized, 
		    "fastafile"     => $fastafile_rna,
		    "dupfile"       => undef,
		    "ebins"         => $ebins,
		    "divbins"       => $divbins,
		    "outfile"       => $expfile_nodups_rna);

      my $cmd = &get_cmd_removedups(\%PARAMS); 
      $pbs->addCmd($cmd);
    }

    #
    # quanrized
    # 
    if ($quantized == 0) {
      $pbs->addCmd("echo \"RNA, quantizing.\"");    
      my $cmd = "perl $scriptdir/quantize_expression_vector.pl -expfile $expfile_nodups_rna -outfile $expfile_q_rna ";	
      if (defined($ebins)) {
	$todo .= " -ebins $ebins ";
      }
      if (defined($divbins)) {
	$todo .= " -divbins $divbins ";
      }
      
      $pbs->addCmd($cmd);
    }


    if ($domifind == 1) { 
      $pbs->addCmd("echo \"RNA, Step 1: seed discovery.\""); 

      my %PARAMS = ("expfile"     => $expfile_nodups_rna,
		    "quantized"   => $quantized,
		    "fastafile"   => $fastafile_rna, 
		    "shuffle"     => $shuffle_mifind, 
		    "kmerfile"    => $kmerfile,
		    "fast"        => $fast,
		    "rna"         => 1,
		    "gap"         => $gap,
		    "k"           => $k,
		    "ebins"       => $ebins,
		    "divbins"     => $divbins,
		    "seedfile"    => $seedfile_rna);
      my $cmd = &get_cmd_mifind(\%PARAMS);
      $pbs->addCmd($cmd);
    }
  
    # if no seeds, exit here
    if ($doskipdiscovery == 0) {
      $pbs->addCmd("CNTM=`perl $scriptdir/count_lines.pl $seedfile_rna`; if [ \$CNTM = 0 ]; then echo \"No DNA motifs. Exiting now.\"; exit; fi");
    }

    if ($domioptimize == 1) {
      $pbs->addCmd("echo \"RNA, Step 2: seed optimization.\""); 
    
      my %PARAMS = ("expfile"     => $expfile_nodups_rna, 
		    "quantized"   => $quantized, 
		    "fastafile"   => $fastafile_rna, 
		    "rna"         => 1, 
		    "gap"         => $gap,
		    "ebins"       => $ebins,
		    "divbins"     => $divbins,
		    "seedfile"    => $seedfile_rna, 
		    "shuffle"     => $shuffle, 
		    "minr"        => $minr, 
		    "add"         => $add,
		    "maxfreq"     => $maxfreq, 
		    "optimslow"   => $optimslow, 
		    "maxdegeneracy" => $maxdegeneracy,		    
		    "optimfile"   => $optimfile_rna);
      my $cmd = &get_cmd_mioptimize(\%PARAMS);
      $pbs->addCmd($cmd);
    }


    if ($doskipdiscovery == 1) {
      $pbs->addCmd("echo \"RNA, Step 2: skip seed optimization, use -motiffile_rna option instead.\"");
      die "Please provide -motiffile_rna\n" if (! -e $motiffile_rna);
      $pbs->addCmd("perl $scriptdir/copy_seeds_to_optim.pl $motiffile_rna $optimfile_rna");
    }

    if ($doskipoptimization == 1) {
      $pbs->addCmd("echo \"RNA, Step 2: skip seed optimization, use seed file instead.\"");
      $pbs->addCmd("perl $scriptdir/copy_seeds_to_optim.pl $seedfile_rna $optimfile_rna");
    }

    

    if ($domisignif == 1) {
      $pbs->addCmd("echo \"RNA, Step 3: evaluation of motif significances.\""); 
    
      my %PARAMS = ("expfile"     => $expfile_nodups_rna, 
		    "quantized"   => $quantized, 
		    "fastafile"   => $fastafile_rna, 
		    "rna"         => 1,
		    "ebins"       => $ebins,
		    "divbins"     => $divbins,
		    "optimfile"   => $optimfile_rna, 
		    "shuffle"     => $shuffle, 
		    "jn"          => $jn, 
		    "jn_t"        => 0, 
		    "jn_f",       => $jn_f,
		    "signiffile"  => $signiffile_rna);
      my $cmd = &get_cmd_misignif(\%PARAMS);
      $pbs->addCmd($cmd);
    }
    
    if ($dobinary == 1) {
      
      my $cmd = "perl $scriptdir/mi_create_overep_expfile.pl -repfile $motifrepfile_rna -expfile $quantized_expfile_rna -quantized $quantized -outexpfile $expfile_b_rna";
      $pbs->addCmd("echo \"RNA, Step 3.5: creating binary expression profiles.\"");
      $pbs->addCmd($cmd);
    }
    
    
  
    if ($domidist == 1) {
      $pbs->addCmd("echo \"RNA, Step 4: discovery of distance constraints.\""); 
    
      my %PARAMS = ("expfile"        => $expfile_nodups_rna, 
		    "expfile_b"      => $expfile_b_rna, 
		    "quantized"      => $quantized, 
		    "fastafile"      => $fastafile_rna, 
		    "rna"            => 1,		  
		    "optimfile"      => $optimfile_rna, 
		    "shuffle"        => $shuffle, 
		    "mbins_interval" => $mbins_interval, 
		    "donbcopies"     => $donbcopies, 
		    "ebins"          => $ebins,
		    "divbins"        => $divbins,
		    "distfile"       => $distfile_rna);
      my $cmd = &get_cmd_midist(\%PARAMS);
      $pbs->addCmd($cmd);
    }
    
    if (($docons == 1) && (defined($fc_fasta1_rna))) {

      $pbs->addCmd("echo \"RNA, Step 5.5: conservation.\"");    
	    
      my $cmd = "perl $scriptdir/evaluate_motifs_conservation.pl -summaryfile $optimfile_rna -fasta1 $fc_fasta1_rna -fasta2 $fc_fasta2_rna -nbgenes $fc_nbgenes_rna -kmerfile $kmerfile_rna -outfile $consfile_rna -rna 1 ";
      if ($shuffle_cons == 1) {
	$cmd .= " -shuffle 1";
      }
      $pbs->addCmd($cmd);
    }
    
    if ($domisummarize == 1) {
      $pbs->addCmd("echo \"RNA, Step 5: summarize information.\"");

      my %PARAMS = ("expfile_rna"        => $expfile_nodups_rna,
		    "fastafile_rna"      => $fastafile_rna,
		    "optimfile_rna"      => $optimfile_rna,
		    "signiffile_rna"     => $signiffile_rna,
		    "distfile_rna"       => $distfile_rna,
		    "consfile_rna"       => (defined($fc_fasta1_rna)?$consfile_rna:undef),
		    "rna"                => 1,
		    "quantized"          => $quantized,
		    "minrobustness"      => $jn_t,
                    "maxrobustness"      => $jn,
		    "removecols_draw"    => $removecols_draw,
		    "sortrowsbyphase"    => $sortmotifsbyphase,
		    "summaryfile"        => $summaryfile_rna,
		    "matrixfile"         => $matrixfile_rna,
		    "fullmatrixfile"     => $fullmatrixfile_rna,
		    "densityfile"        => $densityfile_rna,
		    "columnsfile"        => $columnsfile_rna,
		    "signifcolumnsfile"  => $signifcolumnsfile_rna,
		    "outlevel"           => $outlevel,
		    "minp"               => $minp,
		    "minz"               => $minz,	
		    "oribiasonly"        => $oribiasonly,
		    "namesfile"          => $namesfile_rna,
		    
		    "expfile_b_rna"      => $expfile_b_rna,
		    "outexpfile_b"       => $expfile_b_summary_rna
		   );

      my $cmd = &get_cmd_misummarize(\%PARAMS);
      $pbs->addCmd($cmd);
  
    }

    
    if ($dogetprofiles == 1) {
      
      $pbs->addCmd("echo \"RNA, Step 5.3: generate profiles for motifs in summary file.\"");
      my $cmd = "perl $scriptdir/generate_motif_profiles.pl -fastafile $fastafile_rna -summaryfile $summaryfile_rna -outfile $profiles_rna -rna 1";
      $pbs->addCmd($cmd);

       if (defined($fc_fasta2_rna)) {
	$pbs->addCmd("echo \"RNA, Step 5.3 bis: Generate profiles for $fc_fasta2_rna.\"");
	my $cmd = "perl $scriptdir/generate_motif_profiles.pl -fastafile $fc_fasta2_rna -summaryfile $summaryfile_rna -outfile $orth_profiles_rna -rna 1";
	$pbs->addCmd($cmd);
      }

    }


    if ($dotargets == 1) {
      
      $pbs->addCmd("echo \"RNA, Determining target gene sets.\"");
      my $cmd = "perl $scriptdir/get_motif_targets.pl --matrixfile=$fullmatrixfile_rna --profiles=$profiles_rna --expfile=$quantized_expfile_rna --outfile=$targets_rna";
      $pbs->addCmd($cmd);

    }

    
    if ($dohistograms == 1) {

      $pbs->addCmd("echo \"RNA, Drawing position histograms.\"");
      my $cmd = " perl $scriptdir/draw_position_histogram.pl  --profiles=$profiles_rna --expfile=$quantized_expfile_rna --outdir=$histogramdir_rna --targetfile=$targets_rna --rna=0 --seqlen=$seqlen_rna --bins=10 ";
      $pbs->addCmd($cmd);

    }

    
  
    if (($dornacomp == 1) && ($domicombine == 1)) {

      $pbs->addCmd("echo \"RNA, Step 6: cluster motifs.\"");
    
      my %PARAMS = ("expfile"               => $expfile_b_summary_rna, 
		    "quantized"             => 1,
		    "fastafile_rna"         => $fastafile_rna, 
		    "summaryfile"           => $summaryfile_rna,
		    "poscor"                => $poscor, 
		    "doallstats"            => $doallstats,
		    "clusterfile"           => $clusterfile_rna,
		    "mimatrixfile"          => $mimatrixfile_rna, 
		    "fullmimatrixfile"      => $fullmimatrixfile_rna);

      my $cmd = &get_cmd_micombine(\%PARAMS);
      $pbs->addCmd($cmd);

         
      #
      # do GO analysis of the overlap set between motifs 
      #
      #if (($dogomotifs == 1) && ($dodnarna == 0)) {
#	$cmd = "perl  $scriptdir/motif_pair_go_enrichment.pl --expfile=$expfile_nodups_rna --micombinefile=$mimatrixfile_rna --profiles_dna=$profiles_rna --goindex=$goindex --gonames=$gonames --outfile=$gofile_pairs_rna --outgo=$gofile_pairs_full_rna --maxcatsize=$maxgocatsize ";
	
#	#
#	#  fafner code to fix malfunctioning Hypergeom::cumhyper on fafner only
#	#
#	if (defined($platform) && ($platform eq "fafner")) {
#	  $cmd .= " --usemodule=0 ";
#	}
	
#	$pbs->addCmd($cmd);
#      }

    }

   

    if ($dogomotifs == 1) {
      
      $pbs->addCmd("echo \"RNA, Step 5.6: generate genome-wide GO enrichments for motifs in summary file.\"");

      my $cmd = "perl $scriptdir/motif_single_go_enrichment.pl --summaryfile=$summaryfile_rna --profiles=$profiles_rna --goindex=$goindex --gonames=$gonames --outgo=$gomofile_rna --outfile=$gomofile_full_rna --expfile=$quantized_expfile_rna --maxcatsize=$maxgocatsize ";

      #
      #  fafner code to fix malfunctioning Hypergeom::cumhyper on fafner only
      #
      if (defined($platform) && ($platform eq "fafner")) {
	$cmd .= " --usemodule=0 ";
      }

      if (($dornacomp == 1) && ($domicombine == 1)) {
	$cmd .= " --clusterfile=$clusterfile_rna";
      }
      
      $pbs->addCmd($cmd);

      
    }


    if ($dodrawmotifmaps == 1) {
      #
      # draw interaction motif maps
      #
      $pbs->addCmd("echo \"RNA, draw motif maps for co-localizing pairs\"");
      mkdir "$mimatrixfile_rna\_OUT" if (! -e "$mimatrixfile_rna\_OUT");

      my $cmd = "perl $scriptdir/mi_draw_all_pair_motif_maps.pl --expfile=$quantized_expfile_rna --profiles=$profiles_rna --seqlen=$seqlen_rna --rna=1 --summaryfile=$summaryfile_rna --leftlabel=$leftlabel_rna --rightlabel=$rightlabel_rna --rootdir=$scriptdir --fullmatrixfile=$fullmatrixfile_rna ";

      if (($dornacomp == 1) && ($domicombine == 1)) {
        $cmd .= " --mimatrixfile=$mimatrixfile_rna ";
      }

      $pbs->addCmd($cmd);
    }

    


    if (($dogoclusters == 1) && (defined($goindex))) {
      
      $pbs->addCmd("echo \"RNA, Step 6.5: GO analysis.\"");    

      my $cmd = "perl $scriptdir/clusters_go_enrichment.pl --clusters=$quantized_expfile_rna --goindex=$goindex --gonames=$gonames --N=-1 --outfile=$gofile_full_rna --outgo=$gofile_rna --maxcatsize=$maxgocatsize";

      #
      #  fafner code to fix malfunctioning Hypergeom::cumhyper on fafner only
      #
      if (defined($platform) && ($platform eq "fafner")) {
	$cmd .= " --usemodule=0";
      }

      $pbs->addCmd($cmd);    

    }

       
    if ($domotifnames == 1) {
      
      if (defined($acelibrary_rna)) {
	my $namesfile_rna_tmp = "$namesfile_rna.tmp";
	$pbs->addCmd("echo \"RNA, Step 6.7: identifying motifs.\"");
	my $cmd = "perl $scriptdir/annotate_motifs_using_ace_library.pl --acefiles=\"$acelibrary_rna\" --summaryfile=$summaryfile_rna --namefile=$namesfile_rna --rna=1 --outfile=$namesfile_rna_tmp";	
	$pbs->addCmd($cmd);
	$pbs->addCmd("cp $namesfile_rna $namesfile_rna.old");
	$pbs->addCmd("cp $namesfile_rna.tmp $namesfile_rna");
      }

      if (defined($micrornas)) {
	my $namesfile_rna_tmp = "$namesfile_rna.tmp";
	$pbs->addCmd("echo \"RNA, Step 6.8: identifying motifs that match micrornas.\"");
	my $cmd = "perl $scriptdir/annotate_motifs_using_micrornas.pl --micrornas=$micrornas --summaryfile=$summaryfile_rna --namefile=$namesfile_rna > $namesfile_rna_tmp";	
	$pbs->addCmd($cmd);
	$pbs->addCmd("cp $namesfile_rna $namesfile_rna.old_bef_mirnas");
	$pbs->addCmd("cp $namesfile_rna.tmp $namesfile_rna");
      }
      
    }


    if (($dornacomp == 1) && ($domidrawmatrix == 1)) {

      $pbs->addCmd("echo \"RNA, Step 7: draw matrix figure.\"");

      # columnsfile used only for GO analysis
      my %PARAMS = ("matrixfile"      => $matrixfile_rna,
		    "summaryfile"     => $summaryfile_rna,
		    "columnsfile"     => $columnsfile_rna,
		    "expfile"         => $quantized_expfile_rna,
		    "gofile"          => (defined($goindex)?$gofile_rna:undef),
		    "ps2pdf"          => 1,
		    "every"           => 1,
		    "lp_t_draw"       => $lp_t_draw,
		    "quantized"       => $quantized,
		    "colmap"          => "$scriptdir/HEATMAPS/cmap2.txt",
		    "namesfile"       => $namesfile_rna);
      $PARAMS{"clusterfile"} = $clusterfile_rna if (($sortmotifsbyphase == 0) && ($domicombine >= 1));      

      my $cmd  = &get_cmd_midrawmatrix(\%PARAMS);
      $pbs->addCmd($cmd);


      $pbs->addCmd("echo \"RNA, Step 7b: draw density figure.\"");

      # matrix
      my %PARAMS = ("matrixfile"      => $matrixfile_rna,
                    "densityfile"     => $densityfile_rna,		    
		    "summaryfile"     => $summaryfile_rna,		    
		    "columnsfile"     => $columnsfile_rna,
		    "expfile"         => $quantized_expfile_rna,
		    "gofile"          => (defined($goindex)?$gofile_rna:undef),
		    "ps2pdf"          => 1,
		    "every"           => 1,
		    "quantized"       => $quantized,
		    "redoweblogo"     => 0,
		    "colmap"          => "$scriptdir/HEATMAPS/cmap_dens.txt",
		    "outeps"          => $epsdensityfile_rna,		    
		    "namesfile"       => $namesfile_rna);
      

      $PARAMS{"clusterfile"} = $clusterfile_rna if (($sortmotifsbyphase == 0) && ($domicombine >= 1));      

      my $cmd  = &get_cmd_midrawmatrix(\%PARAMS);
      $pbs->addCmd($cmd);

    }

  
    if ($domidrawdist == 1) {
  
      my $todo = "perl  $scriptdir/mi_create_dist_pvalues.pl -expfile $quantized_expfile_rna -profiles $profiles_rna -distrepfile $distfile_rna.rep -summaryfile $summaryfile_rna -motifdir $summaryfile_rna\_OUT -quantized $quantized -seqlen $seqlen_rna -rightlabel $rightlabel_rna -leftlabel $leftlabel_rna -outdistmatrix $outdistmatrix_rna  -outoriematrix $outoriematrix_rna ";
      $pbs->addCmd($todo);
      
    }


    if (($dornacomp == 1) && ($domidrawinteractions == 1)) {
      my %PARAMS= ("summaryfile"       => $summaryfile_rna,
		   "fullmimatrixfile"  => $fullmimatrixfile_rna,
		   "mimatrixfile"      => $mimatrixfile_rna,
		   "motifdir"          => "$summaryfile_rna\_OUT");
      $PARAMS{"clusterfile"} = $clusterfile_rna if (($sortmotifsbyphase == 0) && ($domicombine >= 1));      
      if ($sortmotifsbyphase == 1) {
	$PARAMS{"orderfile"} = $matrixfile_rna;
      }
      my $cmd = &get_cmd_midrawinteractions(\%PARAMS);
      $pbs->addCmd($cmd);
    
    }

    if ($domotifreport == 1) {
      
      $pbs->addCmd("echo \"RNA, Step 10: motif report.\"");

      my $cmd = "perl $scriptdir/generate_motif_occurrence_summary.pl -expfile $expfile_nodups_rna -summaryfile $summaryfile_rna -profiles $profiles_rna -rna 1  -pvmatrix $matrixfile_rna -columnfile $columnsfile_rna -distmatrix $outdistmatrix_rna -oriematrix $outoriematrix_rna -quantized $quantized ";
      $cmd .= " -gomotifs $gomofile_rna -goindex $goindex " if (defined($goindex));
      $cmd .= " -profiles_orth $orth_profiles_rna " if (defined($fc_fasta2_rna));
      $cmd .= " -mapping_orth $mapping_orth " if (defined($mapping_orth));
      
      if ($quantized == 0) {
	$cmd .= " -expfile_q $expfile_q_rna ";
      }

      $cmd .= " > $motifreport_rna";
      $pbs->addCmd($cmd);
  
    }
    
    
 #   if ($doreduce == 1) {    
#      my $td = "$expfile_nodups_rna.REDUCE_OUT";
#      mkdir $td if (! -e $td);
#      my $cmd = "MATRIXREDUCE/Linux/MatrixREDUCE/bin/MatrixREDUCE -sequence=$fastafile_rna -expression=$expfile_nodups_rna -output=$td -runlog=stdout  -p_value=0.0001 -dyad_length=3 -min_gap=1 -max_gap=1 -flank=1 -single_strand > $expfile_nodups_rna.REDUCE";  
#      $pbs->addCmd($cmd);
#    }
    
    $pbs->addCmd("date");

    
    if ($debug == 1) {
       $pbs->print;
    }

    if ($submit == 0) {
      $pbs->execute;
    } else {
      $jobid_rna = $pbs->submit; print "Submitted job $jobid_dna.\n";
    }
    
  }

  #
  #   end RNA RNA RNA RNA 
  #


  #
  #   combine DNA / RNA
  #

  my $target_dir_dna_rna = "$target_dir/DNA_RNA";
  if (! -e $target_dir_dna_rna) {
    mkdir $target_dir_dna_rna;
  }

  my $expfile_nodups_dna_rna      = "$target_dir_dna_rna/$expfile_file";
  my $expfile_q_dna_rna           = "$expfile_nodups_dna_rna.quantized";
  my $quantized_expfile_dna_rna   = $expfile_nodups_dna_rna;
  if ($quantized == 0) {
    $quantized_expfile_dna_rna = $expfile_q_dna_rna;
  }


  my $summaryfile_dna_rna         = "$expfile_nodups_dna_rna.summary";
  my $densityfile_dna_rna         = "$expfile_nodups_dna_rna.densities";
  my $epsdensityfile_dna_rna      = "$expfile_nodups_dna_rna.densities.eps";
  my $matrixfile_dna_rna          = "$expfile_nodups_dna_rna.matrix";
  my $fullmatrixfile_dna_rna      = "$expfile_nodups_dna_rna.fullmatrix";
  my $columnsfile_dna_rna         = "$expfile_nodups_dna_rna.columns";
  my $signifcolumnsfile_dna_rna   = "$expfile_nodups_dna_rna.signifcolumns";
  my $namesfile_dna_rna           = "$expfile_nodups_dna_rna.motifnames";
  my $clusterfile_dna_rna         = "$expfile_nodups_dna_rna.clusters";
  my $mimatrixfile_dna_rna        = "$expfile_nodups_dna_rna.mimatrix";
  my $fullmimatrixfile_dna_rna    = "$expfile_nodups_dna_rna.fullmimatrix";
  my $gofile_dna_rna              = "$expfile_nodups_dna_rna.GO";
  my $gofile_full_dna_rna         = "$expfile_nodups_dna_rna.GO.full";
  my $gofile_pairs_dna_rna        = "$expfile_nodups_dna_rna.GOpairs";
  my $gofile_pairs_full_dna_rna   = "$expfile_nodups_dna_rna.GOpairs_full";
  my $gomofile_dna_rna            = "$expfile_nodups_dna_rna.GOmotifs";
  my $new_signiffile_dna          = "$expfile_nodups_dna_rna.signif.dna";
  my $new_signiffile_rna          = "$expfile_nodups_dna_rna.signif.rna";
  my $expfile_b_summary_dna_rna   = "$expfile_nodups_dna_rna.summary.binary";

  if ($dodnarna == 1) {

    #
    # make sure sequence file exists
    #
    die "No DNA or no RNA sequence data. Perhaps you need to download the relevant 
species-specific data file from http://tavazoielab.princeton.edu/FIRE.\n" if ((! -e $fastafile_rna) || (! -e $fastafile_dna)); 
 
    #
    # start script
    #
    my $pbs = PBS->new;

    $pbs->setPlatform($platform) if (defined($platform));
    $pbs->setScriptName("$expfile_nodups_dna_rna.script");
    
    $pbs->setQueue($queue)       if (defined($queue));


    $pbs->addCmd("echo \"Processing $expfile_nodups_dna_rna\"");

    $pbs->setWallTime($walltime);

    $pbs->addCmd("date");
    $pbs->addCmd("cd $pwd");

    if ($platform eq 'tcluster') {
      $pbs->addCmd("setenv FIREDIR $firedir");	      
      $pbs->addCmd("setenv LD_LIBRARY_PATH$ENV{FIREDIR}/modules/lib");      
    } else {
      $pbs->addCmd("export FIREDIR=$firedir");	
      $pbs->addCmd("export LD_LIBRARY_PATH=$ENV{FIREDIR}/modules/lib");
    }

    #if (defined($platform) && ($platform eq "fafner")) {
    #  $pbs->addCmd("export DYLD_LIBRARY_PATH=/Genomics/fafner/grid/users/elemento/usr/lib");
    #}

    
    # if no seeds, exit here    
    if ($doskipdiscovery == 0) {      
      $pbs->addCmd("CNTM=`perl $scriptdir/count_lines.pl $seedfile_dna $seedfile_rna`; if [ \$CNTM = 0 ]; then echo \"No DNA or RNA motifs. Exiting now.\"; exit; fi");
    }
    
    if ($quantized == 0) {

      #
      # must find the intersection expression dataset
      #
      my $cmd = "perl $scriptdir/intersect_expression_vectors.pl -expfile1 $expfile_nodups_dna -expfile2 $expfile_nodups_rna -outfile $expfile_nodups_dna_rna";
      $pbs->addCmd($cmd);

      
      #
      # quantize
      #
      if ($quantized == 0) {
	$pbs->addCmd("echo \"DNA/RNA, Step 6.2: quantizing.\"");    
	my $cmd = "perl $scriptdir/quantize_expression_vector.pl -expfile $expfile_nodups_dna_rna -outfile $expfile_q_dna_rna ";	
	if (defined($ebins)) {
	  $cmd .= " -ebins $ebins ";
	}
	if (defined($divbins)) {
	  $cmd .= " -divbins $divbins ";
	}
	
	$pbs->addCmd($cmd);
      }

      #
      # must redo the mi_signif reports ONLY
      #      
      if ($domisignif == 1) {
	
	$pbs->addCmd("echo \"DNA/RNA, Step 3: re-generating motif reports for DNA.\""); 	
	my %PARAMS = ("expfile"     => $expfile_nodups_dna_rna, 
		      "quantized"   => $quantized, 
		      "fastafile"   => $fastafile_dna, 
		      "rna"         => 1,
		      "ebins"       => $ebins,
		      "divbins"     => $divbins,
		      "optimfile"   => $optimfile_dna, 
		      "shuffle"     => $shuffle, 
		      "jn"          => $jn, 
		      "jn_t"        => 0, 
		      "jn_f",       => $jn_f,
		      "doreportonly"=> 1,
		      "signiffile"  => $new_signiffile_dna);

	

	my $cmd = &get_cmd_misignif(\%PARAMS);
	$pbs->addCmd($cmd);

	$pbs->addCmd("echo \"DNA/RNA, Step 3: re-generating motif reports for RNA.\""); 	
	my %PARAMS = ("expfile"     => $expfile_nodups_dna_rna, 
		      "quantized"   => $quantized, 
		      "fastafile"   => $fastafile_rna, 
		      "rna"         => 1,
		      "ebins"       => $ebins,
		      "divbins"     => $divbins,
		      "optimfile"   => $optimfile_rna, 
		      "shuffle"     => $shuffle, 
		      "jn"          => $jn, 
		      "jn_t"        => 0, 
		      "jn_f",       => $jn_f,
		      "doreportonly"=> 1,
		      "signiffile"  => $new_signiffile_rna);
	
	

	my $cmd = &get_cmd_misignif(\%PARAMS);
	$pbs->addCmd($cmd);
      }
      
      
    }

    
    if ($domisummarize == 1) {
      
      $pbs->addCmd("echo \"RNA, DNA Step 5: summarize information.\"");

      my %PARAMS = (
		    "expfile_dna"        => $expfile_nodups_dna,
		    "fastafile_dna"      => $fastafile_dna,
		    "optimfile_dna"      => $optimfile_dna,
		    "signiffile_dna"     => $signiffile_dna,
		    "motifrepfile_dna"   => ($quantized==0?"$new_signiffile_dna.motifs.rep":undef),
		    "distfile_dna"       => $distfile_dna,
		    "consfile_dna"       => (defined($fc_fasta1_dna)?$consfile_dna:undef),
		
		    "expfile_rna"        => $expfile_nodups_rna,
		    "fastafile_rna"      => $fastafile_rna,
		    "optimfile_rna"      => $optimfile_rna,
		    "signiffile_rna"     => $signiffile_rna,
		    "motifrepfile_rna"   => ($quantized==0?"$new_signiffile_rna.motifs.rep":undef),
		    "distfile_rna"       => $distfile_rna,
		    "consfile_rna"       => (defined($fc_fasta1_rna)?$consfile_rna:undef),
		    "rna"                => 2,
		    "quantized"          => $quantized,
		    "minrobustness"      => $jn_t,
                    "maxrobustness"      => $jn,
		    "removecols_draw"    => $removecols_draw,
		    "sortrowsbyphase"    => $sortmotifsbyphase,
		    "outlevel"           => $outlevel,
		    "minp"               => $minp,
		    "minz"               => $minz,
		    "oribiasonly"        => $oribiasonly,
		    "outexpfile"         => $expfile_nodups_dna_rna,
		    "summaryfile"        => $summaryfile_dna_rna,
		    "matrixfile"         => $matrixfile_dna_rna,
		    "fullmatrixfile"     => $fullmatrixfile_dna_rna,	
		    "densityfile"        => $densityfile_dna_rna,
		    "columnsfile"        => $columnsfile_dna_rna,
		    "signifcolumnsfile"  => $signifcolumnsfile_dna_rna,
		    "namesfile"          => $namesfile_dna_rna,

		    "expfile_b_dna"      => $expfile_b_dna,
		    "expfile_b_rna"      => $expfile_b_rna,
		    "outexpfile_b"       => $expfile_b_summary_dna_rna
		   );

      my $cmd = &get_cmd_misummarize(\%PARAMS);
      $pbs->addCmd($cmd);
  
    }
  
    if (($domicombine == 1) && ($sortmotifsbyphase == 0)) {

      $pbs->addCmd("echo \"RNA, Step 6: cluster motifs.\"");
    
      my %PARAMS = ("expfile"               => $expfile_b_summary_dna_rna,
		    "quantized"             => 1,
		    "fastafile_dna"         => $fastafile_dna,
		    "fastafile_rna"         => $fastafile_rna, 
		    "summaryfile"           => $summaryfile_dna_rna,
		    "poscor"                => $poscor, 
		    "doallstats"            => $doallstats,
		    "clusterfile"           => $clusterfile_dna_rna,
		    "mimatrixfile"          => $mimatrixfile_dna_rna, 
		    "fullmimatrixfile"      => $fullmimatrixfile_dna_rna);

      my $cmd = &get_cmd_micombine(\%PARAMS);
      $pbs->addCmd($cmd);

         
      #
      # do GO analysis of the overlap set between motifs 
      #
      #if ($dogomotifs == 1) { 
#	$cmd = "perl  $scriptdir/motif_pair_go_enrichment.pl --expfile=$expfile_nodups_dna_rna --micombinefile=$mimatrixfile_dna_rna --profiles_dna=$profiles_dna --profiles_rna=$profiles_rna --goindex=$goindex --gonames=$gonames --outfile=$gofile_pairs_dna_rna --outgo=$gofile_pairs_full_dna_rna --maxcatsize=$maxgocatsize ";

#        #
#        #  fafner code to fix malfunctioning Hypergeom::cumhyper on fafner only
#        #
#        if (defined($platform) && ($platform eq "fafner")) {
#	  $cmd .= " --usemodule=0 ";
#        }

#	$pbs->addCmd($cmd);
#      }

    }
    


    if (($dogoclusters == 1) && (defined($goindex))) {

      $pbs->addCmd("echo \"DNA/RNA, Step 6.5: GO analysis.\"");    
      my $cmd = "perl $scriptdir/clusters_go_enrichment.pl --clusters=$quantized_expfile_dna_rna --goindex=$goindex --gonames=$gonames --N=-1 --outfile=$gofile_full_dna_rna --outgo=$gofile_dna_rna --maxcatsize=$maxgocatsize ";

      #
      #  fafner code to fix malfunctioning Hypergeom::cumhyper on fafner only
      #
      if (defined($platform) && ($platform eq "fafner")) {
	$cmd .= " --usemodule=0 ";
      }

      $pbs->addCmd($cmd);    

    }

    
    #
    # create name file from the dna/rna ones
    #
    $pbs->addCmd("perl $scriptdir/cat_ifexist.pl $namesfile_dna $namesfile_rna > $namesfile_dna_rna");
    

    #
    # create GO file from the dna/rna ones
    #
    $pbs->addCmd("perl $scriptdir/cat_ifexist.pl $gomofile_dna $gomofile_rna > $gomofile_dna_rna");
    


    if ($domidrawmatrix == 1) {
      $pbs->addCmd("echo \"DNA,RNA, Step 7: draw matrix figure.\"");

      # columnsfile used only for GO analysis
      my %PARAMS = ("matrixfile"      => $matrixfile_dna_rna,
		    "summaryfile"     => $summaryfile_dna_rna,
		    "columnsfile"     => $columnsfile_dna_rna,
		    "expfile"         => $quantized_expfile_dna_rna,		    
		    "gofile"          => (defined($goindex)?$gofile_dna_rna:undef),
		    "ps2pdf"          => 1,
		    "every"           => 1,
		    "lp_t_draw"       => $lp_t_draw,
		    "quantized"       => $quantized,
		    "colmap"          => "$scriptdir/HEATMAPS/cmap2.txt",
		    "namesfile"       => $namesfile_dna_rna);
      #$PARAMS{"clusterfile"} = $clusterfile_dna_rna if ($sortmotifsbyphase == 0);      
      $PARAMS{"clusterfile"} = $clusterfile_dna_rna if (($sortmotifsbyphase == 0) && ($domicombine >= 1));      

      my $cmd  = &get_cmd_midrawmatrix(\%PARAMS);
      $pbs->addCmd($cmd);


      $pbs->addCmd("echo \"DNA,RNA, Step 7: draw matrix density figure.\"");
      
      # density matrix
      my %PARAMS = ("matrixfile"      => $matrixfile_dna_rna,
                    "densityfile"     => $densityfile_dna_rna,		    
		    "summaryfile"     => $summaryfile_dna_rna,
		    "columnsfile"     => $columnsfile_dna_rna,
		    "expfile"         => $quantized_expfile_dna_rna,
		    "gofile"          => (defined($goindex)?$gofile_dna_rna:undef),
		    "ps2pdf"          => 1,
		    "every"           => 1,
		    "quantized"       => $quantized,
		    "redoweblogo"     => 0,
		    "colmap"          => "$scriptdir/HEATMAPS/cmap_dens.txt",
		    "outeps"          => $epsdensityfile_dna_rna,
		    "namesfile"       => $namesfile_dna_rna);
      #$PARAMS{"clusterfile"} = $clusterfile_dna_rna if ($sortmotifsbyphase == 0); 
     
      $PARAMS{"clusterfile"} = $clusterfile_dna_rna if (($sortmotifsbyphase == 0) && ($domicombine >= 1));      
      


      my $cmd  = &get_cmd_midrawmatrix(\%PARAMS);
      $pbs->addCmd($cmd);


    }


    if ($dodrawmotifmaps == 1) {
      #
      # draw interaction motif maps
      #
      $pbs->addCmd("echo \"DNA/RNA, draw motif maps for single motifs and co-localizing pairs\"");
      mkdir "$mimatrixfile_dna_rna\_OUT" if (! -e "$mimatrixfile_dna_rna\_OUT");

      #
      # do DNA first
      # 
      my $cmd = "perl $scriptdir/mi_draw_all_pair_motif_maps.pl --expfile=$quantized_expfile_dna_rna --profiles=$profiles_dna --seqlen=$seqlen_dna --rna=0 --summaryfile=$summaryfile_dna_rna --leftlabel=$leftlabel_dna --rightlabel=$rightlabel_dna --rootdir=$scriptdir --mimatrixfile=$mimatrixfile_dna_rna --fullmatrixfile=$fullmatrixfile_dna_rna "; # > $mimatrixfile_dna\_OUT/index.html";
      $pbs->addCmd($cmd);

      #
      # then do RNA
      # 
      my $cmd = "perl $scriptdir/mi_draw_all_pair_motif_maps.pl --expfile=$quantized_expfile_dna_rna --profiles=$profiles_rna --seqlen=$seqlen_rna --rna=1 --summaryfile=$summaryfile_dna_rna --leftlabel=$leftlabel_rna --rightlabel=$rightlabel_rna --rootdir=$scriptdir --mimatrixfile=$mimatrixfile_dna_rna --fullmatrixfile=$fullmatrixfile_dna_rna "; # > $mimatrixfile_dna\_OUT/index.html";
      $pbs->addCmd($cmd);

      
    }


  
    if ($domidrawinteractions == 1) {
      my %PARAMS= ("summaryfile"       => $summaryfile_dna_rna,
		   "fullmimatrixfile"  => $fullmimatrixfile_dna_rna,
		   "mimatrixfile"      => $mimatrixfile_dna_rna,
		   "motifdir"          => "$summaryfile_dna_rna\_OUT");
      $PARAMS{"clusterfile"} = $clusterfile_dna_rna if ($sortmotifsbyphase == 0);
      if ($sortmotifsbyphase == 1) {
	$PARAMS{"orderfile"} = $matrixfile_dna_rna;
      }
      my $cmd = &get_cmd_midrawinteractions(\%PARAMS);
      $pbs->addCmd($cmd);
    
    }
    
	$pbs->addCmd("date");

	
    if ($debug == 1) {
       $pbs->print;
    }

    if ($submit == 0) {
      $pbs->execute;
    } else {
      
      if (defined($jobid_dna)) {
	$pbs->addDepJob($jobid_dna);
      }
      if (defined($jobid_rna)) {
	$pbs->addDepJob($jobid_rna);
      }
      
      $jobid_dnarna = $pbs->submit;  print "Submitted job $jobid_dnarna.\n";
      
    }
	      
  }

  #
  #   END combine DNA / RNA
  #

}



sub get_cmd_removedups {

  my ($p) = @_;
  
  my $todo = "perl $scriptdir/remove_homologous_sequences_withseed.pl -expfile $p->{expfile} -quantized $p->{quantized} -fastafile $p->{fastafile} -outfile $p->{outfile}";
  if (defined($p->{dupfile})) {
    $todo .= " -dupfile $p->{dupfile} ";
  }

  if (defined($p->{ebins})) {
    $todo .= " -ebins $p->{ebins} ";
  }

  if (defined($p->{divbins})) {
    $todo .= " -divbins $p->{divbins} ";
  }
  
  
  return $todo;
  
  
}


sub get_cmd_mifind {
  my ($p) = @_;

  my $todo = "$progdir/mi_find -expfile $p->{expfile} -quantized $p->{quantized} -fastafile $p->{fastafile} -outfile $p->{seedfile} -report 1 -docor 1 -outprm 1 -shuffle $p->{shuffle} -jn 0 -jn_t 0 -verbose 1 -rna $p->{rna}";
  
  if (defined($p->{kmerfile})) {
    $todo .= " -kmerfile $p->{kmerfile} ";
  }
  
  if (defined($p->{ebins})) {
    $todo .= " -ebins $p->{ebins} ";
  }
  
  #if (defined($p->{fast}) && ($p->{fast} == 0)) {
  #  $todo .= " -fast 0 ";
  #}
  
  if (defined($p->{divbins})) {
    $todo .= " -divbins $p->{divbins} ";
  }

  if (defined($p->{maxselected})) {
    $todo .= " -max_selected $p->{maxselected} ";
  }

  if (defined($p->{k})) {
    $todo .= " -k $p->{k} ";
  }

  if (defined($p->{gap})) {
    $todo .= " -gap $p->{gap} ";
  }

  return $todo;

}


sub get_cmd_mioptimize {
  
  my ($p) = @_;
  
  my $todo = "";
  
  if ($p->{optimslow} == 0) {

    $todo = "$progdir/mi_optimize -expfile $p->{expfile} -quantized $p->{quantized} -fastafile $p->{fastafile} -outfile $p->{optimfile} -kmerfile $p->{seedfile} -shuffle $p->{shuffle} -report 0 -log 1 -rna $p->{rna} -verbose 0 -mask 0 -minr $p->{minr} -maxfreq $p->{maxfreq}";

  } else {

    $todo = "$progdir/mi_optimize_slow -expfile $p->{expfile} -quantized $p->{quantized} -fastafile $p->{fastafile} -outfile $p->{optimfile} -kmerfile $p->{seedfile} -shuffle $p->{shuffle} -report 0 -log 1 -rna $p->{rna}";

  }

  if (defined($p->{add})) {
    $todo .= " -add $p->{add} ";
  }
  
  if (defined($p->{gap})) {
    $todo .= " -gap $p->{gap} ";
  }
  
  if (defined($p->{ebins})) {
    $todo .= " -ebins $p->{ebins} ";
  }

  if (defined($p->{divbins})) {
    $todo .= " -divbins $p->{divbins} ";
  }
  
  if (defined($p->{maxdegeneracy})) {
    $todo .= " -maxdegeneracy $p->{maxdegeneracy} ";
  }

  return $todo;
}


sub get_cmd_misignif {
  my ($p) = @_;
  
  my $todo = "$progdir/mi_signif -expfile $p->{expfile} -quantized $p->{quantized} -fastafile $p->{fastafile} -motiffile $p->{optimfile} -shuffle $p->{shuffle} -optimout 1 -outfile $p->{signiffile} -rna $p->{rna} -jn $p->{jn} -jn_t $p->{jn_t} -jn_f $p->{jn_f} ";
  
  if (defined($p->{ebins})) {
    $todo .= " -ebins $p->{ebins} ";
  }

  if (defined($p->{divbins})) {
    $todo .= " -divbins $p->{divbins} ";
  }
    
  if (defined($p->{doreportonly})) {
    $todo .= " -doreportonly $p->{doreportonly} ";
  }
  
  return $todo;

}


sub get_cmd_midist {
  # my ($expfile, $quantized, $fastafile, $rna, $optimfile, $shuffle, $mbins_interval, $donbcopies, $outfile) = @_;
  my ($p) = @_;
  
  my $todo = "$progdir/mi_dist -expfile $p->{expfile} -expfile_b $p->{expfile_b} -quantized $p->{quantized} -fastafile $p->{fastafile} -motiffile $p->{optimfile} -shuffle $p->{shuffle} -outfile $p->{distfile} -mbins_interval $p->{mbins_interval} -rna $p->{rna} -donbcopies $p->{donbcopies}";

  if (defined($p->{ebins})) {
    $todo .= " -ebins $p->{ebins} ";
  }

  if (defined($p->{divbins})) {
    $todo .= " -divbins $p->{divbins} ";
  }
  
  return $todo;
}


sub get_cmd_misummarize {
  my ($p) = @_;

  #my ($expfile_nodups_dna, $fastafile_dna, $optimfile_dna, $signiffile_dna, $distfile_dna, $consfile_dna,
  #    $expfile_nodups_rna, $fastafile_rna, $optimfile_rna, $signiffile_rna, $distfile_rna, $consfile_rna,
  #    $rna, $quantized, $minrobustness,
  #    $removecols_draw, $sortrowsbyphase,
  #    $outsummary, $outmatrix, $outcolumns, $outgoodcolumns, $outnames, $outexpfile) = @_;
  
  
  my $todo = "perl $scriptdir/mi_summarize.pl ";

  if (($p->{rna} == 0) || ($p->{rna} == 2)) {
    $todo .= " -optimfile_dna $p->{optimfile_dna} -statfile_dna $p->{signiffile_dna} -distfile_dna $p->{distfile_dna} ";
  }
  
  if (($p->{rna} == 1) || ($p->{rna} == 2)) {
    $todo .= " -optimfile_rna $p->{optimfile_rna} -statfile_rna $p->{signiffile_rna} -distfile_rna $p->{distfile_rna} ";
  }

  if ($p->{rna} == 2) {
    $todo .= " -expfile_dna $p->{expfile_dna} -expfile_rna $p->{expfile_rna} -outexpfile $p->{outexpfile} ";
  }
  
  $todo .= " -quantized $p->{quantized} -outmatrix $p->{matrixfile} -outfullmatrix $p->{fullmatrixfile} -outdensitymatrix $p->{densityfile} -outcolumns $p->{columnsfile} -outgoodcolumns $p->{signifcolumnsfile} -outnames $p->{namesfile} -outsummary $p->{summaryfile} -removecols_exp 0 -minrobustness $p->{minrobustness} -outlevel $p->{outlevel} -maxrobustness $p->{maxrobustness} ";

  if (defined($p->{removecols_draw})) {
    $todo .= " -removecols_draw $p->{removecols_draw} ";
  } 
    
  if (defined($p->{consfile_dna})) {
    $todo .= " -consfile_dna $p->{consfile_dna} ";
  }
  
  if (defined($p->{consfile_rna})) {
    $todo .= " -consfile_rna $p->{consfile_rna} ";
  }  

  if (defined($p->{motifrepfile_dna})) {
    $todo .= " -motifrepfile_dna $p->{motifrepfile_dna} ";
  }
  
  if (defined($p->{motifrepfile_rna})) {
    $todo .= " -motifrepfile_rna $p->{motifrepfile_rna} ";
  }  
  
  if (defined($p->{minp})) {
    $todo .= " -minp $p->{minp} ";
  }  

  if (defined($p->{minz})) {
    $todo .= " -minz $p->{minz} ";
  }  

  if (defined($p->{oribiasonly})) {
    $todo .= " -oribiasonly $p->{oribiasonly} ";
  }  

  if (defined($p->{expfile_b_dna})) {
    $todo .= " -expfile_b_dna $p->{expfile_b_dna} ";
  }
  
  if (defined($p->{expfile_b_rna})) {
    $todo .= " -expfile_b_rna $p->{expfile_b_rna} ";
  }
  
  if (defined($p->{outexpfile_b})) {
    $todo .= " -outexpfile_b $p->{outexpfile_b} "; 
  }
    
  $todo .= " -sortrowsbyphase $p->{sortrowsbyphase} ";
 
  return $todo;

}


sub get_cmd_micombine {

  #my ($expfile_nodups, $quantized,
  #    $fastafile_dna, $fastafile_rna, $summaryfile,
  #    $poscor, $doallstats,		 
  #    $clusterfile, $mimatrixfile, $fullmimatrix) = @_;

  my ($p) = @_;

  my $todo = "$progdir/mi_combine -expfile $p->{expfile}";
  
  if (defined($p->{fastafile_dna})) {
    $todo .= " -fastafile_dna $p->{fastafile_dna} ";
  }
  
  if (defined($p->{fastafile_rna})) {
    $todo .= " -fastafile_rna $p->{fastafile_rna} ";
  }
  
  $todo .= " -summaryfile $p->{summaryfile} -quantized $p->{quantized} -outfile $p->{clusterfile} -outmimatrix $p->{mimatrixfile} -outfullmimatrix $p->{fullmimatrixfile} -poscor $p->{poscor} -doallstats $p->{doallstats} ";
  
  return $todo;
  
}


sub get_cmd_midrawmatrix {
  my ($p) = @_;

  my $todo = "perl $scriptdir/mi_draw_matrix.pl --expfile=$p->{expfile} --matfile=$p->{matrixfile} --summaryfile=$p->{summaryfile} --columnsfile=$p->{columnsfile} --ps2pdf=$p->{ps2pdf} --every=$p->{every} --quantized=$p->{quantized} --motifnames=$p->{namesfile} --ybase=250";
  
  if (defined($p->{colmap})) {
    $todo .= " --colmap=$p->{colmap} ";
  }

  if (defined($p->{clusterfile})) {
    $todo .= " --clustfile=$p->{clusterfile} ";
  }

  if (defined($p->{densityfile})) {
    $todo .= " --densityfile=$p->{densityfile} ";
  }  

  if (defined($p->{redoweblogo})) {
    $todo .= " --redoweblogo=$p->{redoweblogo} ";
  }  

  if (defined($p->{outeps})) {
    $todo .= " --outeps=$p->{outeps} ";
  }  

  if (defined($p->{gofile})) {
    $todo .= " --gofile=$p->{gofile} ";
  }  

  if (defined($p->{lp_t_draw})) {
    $todo .= " --lp_t_draw=$p->{lp_t_draw} ";
  }

  return $todo;
}


sub get_cmd_midrawinteractions {
  my ($p) = @_;

  my $todo = "perl $scriptdir/mi_draw_interaction_matrix.pl -summaryfile $p->{summaryfile} -matrixfile $p->{fullmimatrixfile} -resmatrixfile $p->{mimatrixfile} -motifdir $p->{motifdir} ";
  if (defined($p->{clusterfile})) {
    $todo .= " -clustfile $p->{clusterfile} ";
  }
  if (defined($p->{orderfile})) {
    $todo .= " -orderfile $p->{orderfile} ";
  }
  return $todo;
}


sub check_input_file {
  my ($expfile, $exptype) = @_;
  
  my $ta = Table->new;
  $ta->loadFile($expfile);
  my $a_ref = $ta->getArray();
  
  my $r = shift @$a_ref; 
  if ($r->[1] =~ /^\d/) {
    print "WARNING: your file might not contain a header line ($r->[1]).\n";
  }
  
  my %H = ();
  my %V = ();
  foreach my $r (@$a_ref) {
    if (defined($H{$r->[0]})) {
      print "Your files contains multiple rows with the same gene id. Please correct that before applying FIRE.\n";
      return 0;
    }
    $H{$r->[0]} ++;
    $V{$r->[1]} = 1;
  } 
  
  
  my @v = values ( %V );
  @v = sort { $a <=> $b } @v;
  my $max = $#v;
  if (($exptype eq 'discrete') && (scalar(@v) != $max+1)) {
    my $n1 = scalar(@v);
    my $n2 = $max + 1;
    die "Problem. Your discrete vector is missing some symbols ($n1 != $n2).\n";
    return 0;
  }
  
  return 1;
}


sub readSpeciesData {
  my ($species) = @_;
  
  my %H = ();
  open IN, "$ENV{FIREDIR}/FIRE_DATA/SPECIES_DATA/$species" or die "No data file for $species.\n";
  while (my $l = <IN>) {
    chomp $l;
    my @a = split /\t/, $l, -1;    
    if ($a[1] =~ /^FIRE_DATA/) {
      $a[1] = "$ENV{FIREDIR}/$a[1]";
    }	
    $H{$a[0]} = $a[1];    
  }  
  close IN;
  
  return \%H;
}
