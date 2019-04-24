#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

my $scriptdir = "$home/PERL_MODULES/SCRIPTS";

use Sets;
use Getopt::Long;
use strict;

my $gsematrix    = undef;
my $platformid   = undef;
my $isclean      = undef;
my $genelist     = undef;
my $platformfile = undef;
my $logtransform = 0;
my $qnorm        = 0;
my $varnorm      = 0;
my $genenames    = undef;
my $negto1       = 0;

if (@ARGV == 0) {
  die "Usage: perl GEO_process_series_matrix.pl --gsematrix=FILE --platformid=STR [ --negto1=INT --platformfile=FILE --isclean=INT --logtransform=1 --varnorm=1 --qnorm=1 ]\n";
}

GetOptions(
	   'platformid=s'   =>  \$platformid,
	   'platformfile=s' =>  \$platformfile,  
	   'gsematrix=s'    =>  \$gsematrix,
	   'isclean=s'      =>  \$isclean,
	   'genelist=s'     =>  \$genelist,
           'negto1=s'       =>  \$negto1,
	   'qnorm=s'        =>  \$qnorm,
	   'varnorm=s'      =>  \$varnorm,
	   'logtransform=s' =>  \$logtransform,
	   'genenames=s'    =>  \$genelist
	  );

my $dirgenes = "/home/elemento/PROGRAMS/PLASMODIUM/MICROARRAYS";

my %GENELISTS = (
		 "py"    => "$dirgenes/py_genes.txt",
		 "pf"    => "$dirgenes/pf_genes.txt",
                 "pf6.0" => "/Users/olivier/PROGRAMS/PLASMODIUM/GENOMES/annotations/pf_annotation_6.0.txt.genenames",
                 "py6.0" =>     "$ENV{HOME}/PROGRAMS/PLASMODIUM/GENOMES/annotations/py_annotation_6.0.txt.genenames",
		 "pb"    => "$dirgenes/pb_genes.txt",
		 "human" => "$ENV{HOME}/PROGRAMS/ChIPseeqer/DATA/hg18/refGene.txt.07Jun2010.NM",
		 "mouse" => "$ENV{HOME}/PROGRAMS/FIRE/FIRE_DATA/MOUSE/SEQUENCES/mouse_u_1000_0.fa.genenames",
		 "human_orf" => "/Users/olivier/PROGRAMS/CANCERGENES/SANGER/refLink.simplified.human.GENES"
		);

my $gseid = $gsematrix;
my $plaid = undef;
my $todo  = undef;

if (!defined($platformfile)) {

  if (!defined($platformid)) {
    
    #
    # determine platform
    #
    $todo = "perl $scriptdir/GEO_get_platform_from_series_matrix.pl $gseid";
    print "$todo\n";
    $plaid = `$todo`;
    $plaid =~ s/\n//g;
    
  } elsif (defined($platformid)) {
    $plaid = $platformid;
  }


  print "\nPlatform ID is $plaid.\n";

  #
  # download platform if needed 
  #

  

  if (! -e $platformfile) {
    print "No $platformfile available, downloading ...\n";
    system("getgeoplatform $plaid");  
  }

  $platformfile = "$plaid.txt";

}


#
# parse platform 
#

open IN, "$platformfile" or die "Cannot open platform file\n";
my $cntl = 0;
while (my $l = <IN>) {
  #$l =~ s/\r//g;
  my @a = split /\t/, $l, -1;
  if ($l !~ /^[\!\#\^]/) {
    print "$l";
    $cntl ++;
    last if ($cntl == 5);
  }
}
close IN;


print "\nPlease help me ... which column contains the gene identifiers ?\n";

my $c = undef;
while (1) {

  $c = <STDIN>; #1; #<STDIN>;
  chomp $c;
  
  print "\nMean this one ?\n";
  
  open IN, $platformfile;
  $cntl = 0;
  while (my $l = <IN>) {
    chomp $l;
    my @a = split /\t/, $l, -1;
    if (($l !~ /^[\!\#\^]/) && ($a[$c] ne "")) {
      print "$a[$c]\n";
      $cntl ++;
      last if ($cntl == 5);
    }
  }
  close IN;


  my $yes = <STDIN>; #""; #<STDIN>;
  chomp $yes;
  if ($yes eq "yes") {
    last;
  }
}
  

print "Making a simplified version of the platform file ... \n";

$todo = "perl $scriptdir/columns.pl 0 $c < $platformfile  > $platformfile.clean";
#print "$todo\n";

system($todo);

if (-e "$platformfile.clean") {
  print "$platformfile.clean created.\n"
}


my $species = undef;

if (!defined($genelist)) { 
  
  print "\nCleaning up platform file using genes for which we have information (sequence, etc) ... \n";
  print "Which species is this ? [". join("/", keys(%GENELISTS)) . "]\n";
  $species = <STDIN>; #'yoelli'; #<STDIN>; 
  chomp $species;


  if (! -e $GENELISTS{$species}) {
    die "No gene list $GENELISTS{$species} for $species.\n";
  } else {
    $genelist = $GENELISTS{$species};
  }

}

$todo = "perl $scriptdir/ids_consolidate_mapping_using_gene_list.pl $platformfile.clean $genelist > $platformfile.clean.consolidated";
#print "$todo\n";
system($todo);
print "Done.\n";


if (!defined($isclean)) {
  print "\nClean up series matrix ... ";
  $todo = "perl $scriptdir/GEO_series_to_matrix.pl $gseid > $gseid.clean";
  system($todo);
  print "Done ($gseid.clean created).\n";
} else {
  system("cp $gseid $gseid.clean") == 0 or die "Could not copy $gseid to $gseid.clean\n";
} 

print "\nTranslating gene ids inside series matrix ... ";
$todo = "perl $scriptdir/translate_column_using_table_column.pl --table=$gseid.clean --col=0 --dict=$platformfile.clean.consolidated --k=0 --v=1 > $gseid.ids";
system($todo);
print "Done.\n";



my $infile  = undef;
my $outfile = "$gseid.rowavg";

print "\nAveraging rows with same id ... ";
$todo = "perl $scriptdir/expression_average_rows_with_same_id.pl $gseid.ids > $outfile";
system($todo);
print "Done.\n";

if ($logtransform == 1) {

  $infile   = $outfile;
  $outfile .= ".log";

  print "\nLog-transforming the expression values ... ";
  $todo = "perl $scriptdir/expression_log_transform_matrix.pl $infile ";
  if ($negto1 == 1) {
    $todo .= " 1 ";
  }
  $todo .= " > $outfile ";
  system($todo);
  print "Done.\n";
}

if ($qnorm == 1) {

  $infile   = $outfile;
  $outfile .= ".qnorm";

  print "\nRunning a quantile normalization using R ... ";

  my $txt = "f <- \"$infile\"
m <- read.csv(f, sep=\"\\t\", row.names=1, header=T, check.names=F)
library(limma)
mn <- normalizeQuantiles(m)
write.table(format(mn, digits=4, trim=T), file=\"$outfile\", sep=\"\\t\", quote=F, col.names=NA, row.names=T)
";
  
  my $tmpfile = Sets::getTempFile("/tmp/Rscript");
  open OUT, ">$tmpfile" or die "Cannot open $tmpfile for writing.\n";
  print OUT $txt;
  close OUT;
  system("R CMD BATCH $tmpfile");

  print "Done.\n";
  
}

if ($varnorm == 1) {

  $infile   = $outfile;
  $outfile .= ".varnorm";

  print "\nVariance normalize matrix ... ";
  $todo = "perl $scriptdir/expression_normalize_rows.pl < $infile > $outfile";
  system($todo);
  print "Done.\n"; 
}


exit;

print "\nCluster ... \n";
$todo = "getnbclusters $gseid.rowavg";
my $nbc = `$todo`; chomp $nbc;




$todo = "/home/elemento/PROGRAMS/CLUSTER/cluster-1.35/src/cluster.exe -cg a -ng -g 2 -k $nbc -r 10 -f $gseid.rowavg";
system($todo) == 0 or die "Could not cluster ... \n";

print "Done (clustered into $nbc clustered)\n";
