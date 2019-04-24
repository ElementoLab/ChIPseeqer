#!/usr/bin/perl

#
# fire_gapped.pl: wrapper script for running and merging gapped FIRE runs. 
# this script was created by Hani Goodarzi.
#
#

use lib "$ENV{FIREDIR}/SCRIPTS";

use strict ;
use Sets ;
use Getopt::Long;
use PBS;

if ((!$ENV{FIREDIR}) || ($ENV{FIREDIR} eq '')) 
{
    print "The FIREDIR environment variable is not set. Please set it using export or setenv (see FIRE tutorial online).\n";
    exit;
}

if (@ARGV == 0) {
    die "Usage: perl fire_gapped.pl --expfile=FILE --exptype=TXT --fastafile_DNA=FILE --fastafile_RNA=FILE --k=6 --gap=X-Y --dophase1=1/0 --dophase2=1/0 --submit=1/0 --k=6\n";
}


my $expfile       = undef ;
my $exptype       = undef ;
my $fastafile_dna = undef ;
my $fastafile_rna = undef ;
my $species       = undef ;
my $gap           = 0 ;
my $gap1          = 0 ;
my $gap2          = 0 ;
my $kmer          = 6 ;
my $dophase1      = -1 ;
my $dophase2      = -1 ;
my $submit        = 0 ;
my $quantized     = 1 ;
my $ebins         = undef;
my $divbins       = 50;
my $mbins         = 2;
my $skipfireruns  = 0 ;
my $kungapped     = undef;
my $jn_t_gapped   = undef;
my $minr          = undef;
my $platform      = undef;
my $suffix        = undef;

my @argv_copy = @ARGV;  # make a copy before it is erased

Getopt::Long::Configure("pass_through");

GetOptions ( 'expfile=s'          => \$expfile,
	     'fastafile_DNA=s'    => \$fastafile_dna,
	     'fastafile_RNA=s'    => \$fastafile_rna,
	     'species=s'          => \$species,
	     'gap=s'              => \$gap,
	     'submit=s'           => \$submit,
	     'skipfireruns=s'     => \$skipfireruns,
	     'exptype=s'          => \$quantized,
	     'kungapped=s'        => \$kungapped,
	     'jn_t_gapped=s'      => \$jn_t_gapped,
	     'ebins=s'            => \$ebins,
	     'minr=s'             => \$minr,
	     'platform=s'         => \$platform,
	     'divbins=s'          => \$divbins,
	     'suffix=s'           => \$suffix,
	     'mbins=s'            => \$mbins);

if (defined($species)) {
    #
    #  read species file
    #
    my $species_data = readSpeciesData($species);
  
    if (!defined($fastafile_dna)) {
	$fastafile_dna = $species_data->{"fastafile_dna"};
	print "fastafile_dna is $fastafile_dna\n\n";
    }  
  
    if (!defined($fastafile_rna)) {
	$fastafile_rna = $species_data->{"fastafile_rna"};
	print "fastafile_rna is $fastafile_rna\n\n";
    }
}

if ($quantized eq "discrete") {
  $quantized = 1;

} else {

  $quantized = 0;
  if (! defined $ebins) {
    my $s = `wc $expfile -l`;
    my ($l_num, $dummy) = split(/\s/, $s);
    $ebins = (0.5 + $l_num) / ( $divbins * $mbins );
  }
}

##################Phase 1####################
my $firedir = $ENV{'FIREDIR'};
my $pwd     = `pwd`; $pwd =~ s/\n//;
my $metadir = $expfile."\_META" ;
mkdir $metadir if (! -d $metadir) ;
system("cp $expfile $metadir");

my $fn = Sets::filename($expfile);

my @dep_jobs;
if ($gap =~ /-/)
{
    $gap =~ /(\S+)-(\S+)/ ;
    $gap1 = $1 ;
    $gap2 = $2 ;
}
else
{
    $gap1=$gap;
    $gap2=$gap;
}
foreach my $g ($gap1..$gap2)
{
  my @arg = @argv_copy;
  
  last if ($skipfireruns==1);
  system("cp $expfile $expfile\_g$g");
  my $t_expfile = "$expfile\_g$g";
  my $pbs = PBS->new ;

  $pbs->setPlatform($platform) if (defined($platform));
  if ($platform eq 'tcluster') {
    $pbs->addCmd("setenv FIREDIR $firedir");	      
  } else {
    $pbs->addCmd("export FIREDIR=$firedir") ;
  }
  $pbs->addCmd("cd $pwd") ;
  
  my $script_file = "$metadir/$fn\_$g.script";
  $pbs->setScriptName("$script_file");
  $pbs->addCmd("echo \"Running FIRE with gap $g\"") ;
  for (my $i=0 ; $i<@arg ; $i++) {

    if ($arg[$i] =~ /--expfile/) {
      $arg[$i] = "--expfile=$t_expfile" ;
    }

    if ($arg[$i] =~ /--gap/) {
      $arg[$i] = "--gap=$g" ;
    }

    if ($arg[$i] =~ /--submit/) {
      $arg[$i] = "--submit=0" ;
    }

    if ($arg[$i] =~ /\-\-skipfireruns/) {
      $arg[$i] =~ s/\-\-skipfireruns=\d//;
    }

    if ($arg[$i] =~ /\-\-kungapped\=\d+/) {
      # that argument does not need to be passed to the FIRE runs
      $arg[$i] = '';
    }

    if ($arg[$i] =~ /\-\-jn_t_gapped/) {
      $arg[$i] = '';
    }
  }

  if (defined($jn_t_gapped) && ($g != 0)) {
    push @arg, "--jn_t=$jn_t_gapped";
  }

  # if $g == 0, we might want to use k=7
  if (($g == 0) && defined($kungapped)) {
    print "g==0, kungapped==$kungapped\n";
    push @arg, "--k=$kungapped";
    } elsif (grep /\-\-k\=/, @arg) {
      print "found --k= ??\n";
      print join("/", @arg) . "\n";
    }
    elsif (! grep /\-\-k\=/, @arg) {
      print "Setting --k=6.\n";
      push @arg, "--k=6";  
    }

    my $args = join(" ", @arg);
    $pbs->addCmd("perl $firedir/fire.pl $args --dodnarna=0") ;

    print("Command: perl fire.pl $args\n\n") ;

    my $fire_jobid ;
    if ($submit==0)
    {
	$pbs->execute ;
    }
    elsif ($submit==1)
    {
	$fire_jobid = $pbs->submit ;
	push(@dep_jobs, $fire_jobid);
	print "Submitted job $fire_jobid.\n";
    }
}

#
##################Phase 2##################
#
my $pbs = PBS->new ;

$pbs->setPlatform($platform) if (defined($platform));

if ($platform eq 'tcluster') {
  $pbs->addCmd("setenv FIREDIR $firedir");	      
} else {
  $pbs->addCmd("export FIREDIR=$firedir") ;
}

$pbs->addCmd("cd $pwd") ;
foreach my $g ($gap1..$gap2)
{
    $pbs->addCmd("rm $expfile\_g$g") ;
}
print ("Command: perl $firedir/SCRIPTS/combine_gapped_motifs.pl --metadir=$metadir --expfile=$expfile --gap=$gap\n\n") ;
$pbs->addCmd("perl $firedir/SCRIPTS/combine_gapped_motifs.pl --metadir=$metadir --expfile=$expfile --gap=$gap") ;
$pbs->setScriptName("$metadir/combine_motifs.script");
$pbs->addCmd("echo \"Combining motifs...\n\"") ;
my $cmb_jobid ;
if ($submit==0)
{
    $pbs->execute ;
}
elsif ($submit==1)
{
    if (@dep_jobs>0)
    {
	foreach my $id (@dep_jobs)
	{
	    $pbs->addDepJob($id);
	}
    }
    $cmb_jobid = $pbs->submit ;
    print "Submitted job $cmb_jobid.\n";
}

$pbs = PBS->new ;

$pbs->setPlatform($platform) if (defined($platform));

if ($platform eq 'tcluster') {
  $pbs->addCmd("setenv FIREDIR $firedir");	      
} else {
  $pbs->addCmd("export FIREDIR=$firedir") ;
}
$pbs->addCmd("cd $pwd") ;
my $file = $expfile;
$expfile = "$metadir/$fn" ;
my $todo = "$firedir/PROGRAMS/mi_motif_compare -expfile $expfile -motiffile $metadir/motifs_dna_redundant -quantized $quantized -outfile $metadir/motifs_dna -fastafile $fastafile_dna";
if (defined($minr)) {
  $todo .= " -minr $minr ";
}
$pbs->addCmd($todo);
my $script_file = "$metadir/$fn\_cdna.script";
$pbs->setScriptName("$script_file");
$pbs->addCmd("echo \"Removing redundancy in DNA motifs\n\"") ;
my $cdna_jobid ;
if ($submit==0)
{
    $pbs->execute ;
}
elsif ($submit==1)
{
    $pbs->addDepJob($cmb_jobid);
    $cdna_jobid = $pbs->submit ;
    print "Submitted job $cdna_jobid.\n";
}

$pbs = PBS->new ;

$pbs->setPlatform($platform) if (defined($platform));

if ($platform eq 'tcluster') {
  $pbs->addCmd("setenv FIREDIR $firedir");	      
} else {
  $pbs->addCmd("export FIREDIR=$firedir") ;
}

$pbs->addCmd("cd $pwd") ;
$todo = "$firedir/PROGRAMS/mi_motif_compare -expfile $expfile -motiffile $metadir/motifs_rna_redundant -quantized $quantized -outfile $metadir/motifs_rna -fastafile $fastafile_rna -rna 1 ";
if (defined($minr)) {
  $todo .= " -minr $minr ";
}
$pbs->addCmd($todo);
$script_file = "$metadir/$fn\_crna.script";
$pbs->setScriptName("$script_file");
$pbs->addCmd("echo \"Removing redundancy in RNA motifs\n\"") ;
my $crna_jobid ;
if ($submit==0)
{
    $pbs->execute ;
}
elsif ($submit==1)
{
    $pbs->addDepJob($cmb_jobid);
    $crna_jobid = $pbs->submit ;
    print "Submitted job $crna_jobid.\n";
}

my $args = join(" ", @argv_copy);
$args =~ s/--gap=\d+-\d+\s//;
$args =~ s/--skipfireruns=\d//;
$args =~ s/--domisearch=\d//;
$args =~ s/--submit=\d//;
$args =~ s/--expfile=\S+//;
$args =~ s/--kungapped=\d+//;
$args =~ s/--jn_t_gapped=\d+//;

$pbs = PBS->new ;

$pbs->setPlatform($platform) if (defined($platform));

if ($platform eq 'tcluster') {
  $pbs->addCmd("setenv FIREDIR $firedir");	      
} else {
  $pbs->addCmd("export FIREDIR=$firedir") ;
}

$pbs->addCmd("cd $pwd") ;

my $cmd = "perl $firedir/fire.pl --expfile=$expfile $args --dodna=1 --dorna=1 --dodnarna=0 --doskipdiscovery=1 --motiffile_dna=$metadir/motifs_dna --motiffile_rna=$metadir/motifs_rna --submit=0 --dodrawmotifmaps=0 --doalignace=0 --domidrawdist=0 --docons=1 --dogoclusters=0 --domidrawinteractions=0 --domidrawmatrix=0 ";
if (defined($ebins)) {
  $cmd .= " --ebins=$ebins";
}

$pbs->addCmd($cmd);
$pbs->addCmd("perl $firedir/SCRIPTS/combine_signif_files.pl --expfile=$file --gap1=$gap1 --gap2=$gap2 --motiffile_dna=$metadir/motifs_dna --motiffile_rna=$metadir/motifs_rna");
$pbs->setScriptName("$metadir/make_summary_file.script");
my $sum_jobid ;
if ($submit==0)
{
    $pbs->execute ;
}
elsif ($submit==1)
{
    $pbs->addDepJob($cdna_jobid);
    $pbs->addDepJob($crna_jobid);

    $sum_jobid = $pbs->submit ;
    print "Submitted job $sum_jobid.\n";
}

$pbs = PBS->new ;

$pbs->setPlatform($platform) if (defined($platform));
if ($platform eq 'tcluster') {
  $pbs->addCmd("setenv FIREDIR $firedir");	      
} else {
  $pbs->addCmd("export FIREDIR=$firedir") ;
}
$pbs->addCmd("cd $pwd") ;
$pbs->addCmd("perl $firedir/fire.pl --expfile=$expfile $args --dodna=0 --dorna=0 --domisignif=0 --doskipdiscovery=1 --motiffile_dna=$metadir/motifs_dna --motiffile_rna=$metadir/motifs_rna --submit=0") ;
$pbs->setScriptName("$metadir/fire_nondis.script");
if ($submit==0)
{
    $pbs->execute ;
}
elsif ($submit==1)
{
    $pbs->addDepJob($sum_jobid);

    my $jobid = $pbs->submit ;
    print "Submitted job $jobid.\n";
}

sub readSpeciesData 
{
    my ($species) = @_;
    
    my %H = ();
    open IN, "$ENV{FIREDIR}/FIRE_DATA/SPECIES_DATA/$species" or die "No data file for $species.\n";
    while (my $l = <IN>) 
    {
	chomp $l;
	my @a = split /\t/, $l, -1;    
	if ($a[1] =~ /^FIRE_DATA/) 
	{
	    $a[1] = "$ENV{FIREDIR}/$a[1]";
	}
	$H{$a[0]} = $a[1];    
    }  
    close IN;
    
    return \%H;
}
