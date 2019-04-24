#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";

use Sets;
use PBS;
use Getopt::Long;
use strict;

if ($ENV{CHIPSEEQERDIR} eq "") {
  die "Please set CHIPSEEQERDIR\n"; 
}

my $filepat  = "s_*.txt";
my $submit   = 2;
my $align    = 0;
my $annotate = 1;
my $readlen  = undef;
my $genome   = undef;
my $k        = 30;

if (@ARGV == 0) {
  die "Args: --genome=FASTA --k=INT --submit=INT\n";
}

GetOptions("genome=s"   => \$genome,
	   "k=s"        => \$k,
           "submit=s"   => \$submit);



my $pbs = PBS->new;
$pbs->setPlatform("panda");
#$pbs->setQueue("nomyrinet");
#$pbs->useAllNode(1);
#$pbs->setWallTime("24:00:00");
$pbs->setScriptName("script.readmapjob");
$pbs->setWallTime("24:00:00");
#$pbs->setMemory("8g");
$pbs->addCmd("export CHIPSEEQERDIR=$ENV{CHIPSEEQERDIR}");
#$pbs->addCmd("cd $ENV{PWD}");

$pbs->addCmd("date");

# in any case, copying SAM file locally
$pbs->addCmd("cp -v $genome \$TMPDIR/");
my $genomefile = Sets::filename($genome);


$pbs->addCmd("cd \$TMPDIR");

my $gtpath = "/oelab01_scratch001/ole2001/PROGRAMS/TOOLS/genometools-1.3.4/bin";

$pbs->addCmd("$gtpath/gt suffixerator -dna -pl -tis -suf -lcp -v -parts 4 -db $genomefile -indexname reads");
$pbs->addCmd("$gtpath/gt tallymer mkindex -mersize $k -minocc 1 -indexname tyr-reads -counts -pl -esa reads -scan");
$pbs->addCmd("$gtpath/gt tallymer search -output qseqnum qpos counts sequence -tyr tyr-reads -q $genomefile > readmap.txt");

$pbs->addCmd("cp -v readmap.txt $ENV{PWD}/");
$pbs->addCmd("rm -Rf \$TMPDIR/*");

if ($submit == 1) {
  my $jobid = $pbs->submit;
} elsif ($submit == 0) {
  $pbs->execute;
} else {
  $pbs->print;
}

  


