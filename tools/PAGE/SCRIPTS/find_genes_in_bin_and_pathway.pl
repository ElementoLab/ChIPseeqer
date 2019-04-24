#!/usr/bin/perl

use lib "$ENV{PAGEDIR}/SCRIPTS";

use Getopt::Long;
use Table;
use Sets;
use strict;

my $expfile    = undef;
my $bin        = undef;
my $pathway    = undef;
my $species    = undef;
my $onecatfile = undef;
my $exptype    = undef;

die "Please specify --expfile=FILE --pathway=STR --bin=INT --species=STR (e.g. human_go) \n" if (@ARGV == 0);


GetOptions("expfile=s"  => \$expfile,
	   "pathway=s"  => \$pathway,
	   "bin=s"      => \$bin,
	   "exptype=s"  => \$exptype,
	   #"pathways=s" => \$species,
	   "onecatfile=s" => \$onecatfile,
	   "species=s"  => \$species);



open IN, "$expfile\_PAGE/info.txt" or die "Cannot open $expfile\_PAGE/info.txt\n";
my @lines = <IN>;
close IN;

my @go = grep /goindexfile/, @lines;

if (!defined($species)) {
  ($species) = $go[0] =~ /ANNOTATIONS\/(.+?)\//;
  if ($species eq "") {
    die "Problem detecting the pathway database; use --species=pathways (eg human_go_orf) to correct\n";
  }
}

if ( $species =~ /kegg/i){
	$species= "kegg";
}

#
# load pathways
#
my $speciesdir  = "$ENV{PAGEDIR}/PAGE_DATA/ANNOTATIONS/$species";

# load annotation
my $goindexfile = undef;
my $gonamesfile = undef;
my $descfile    = undef;

if (!defined($onecatfile)) {
  $goindexfile = "$speciesdir/$species\_index.txt";
  $gonamesfile = "$speciesdir/$species\_names.txt";
  $descfile    = "$speciesdir/$species\_genedesc.txt";
} else {
  $goindexfile = $onecatfile;
  $gonamesfile = "$onecatfile.names";
}

#$descfile    = "$ENV{PAGEDIR}/PAGE_DATA/ANNOTATIONS/human_go_orf/human_go_orf_genedesc.txt";

my @a_bins = split /\,/, $bin;

# get matching pathway
my $ta = Table->new;
$ta->loadFile($gonamesfile);
my $a_ref = $ta->getArray();
my %matchinggo = ();
foreach my $r (@$a_ref) {
  if (($r->[1] =~ /$pathway/) || ($r->[0] eq $pathway)) {
    $matchinggo{ $r->[0] } = 1;
  }
}

# get genes in the corresp pathway
$ta->loadFile($goindexfile);
$a_ref = $ta->getArray();

my %genes = ();
foreach my $r (@$a_ref) {
  my $g = shift @$r;
  
  foreach my $c (@$r) {

    if (defined($matchinggo{ $c })) {
      $genes{ $g } = 1; 
    }	

  }
  
}

# load $descfile

my $h_ref_desc = undef;

if (-e $descfile) {
  $ta->loadFile($descfile);
  $h_ref_desc = $ta->getIndex(0);
}


# now go thru the expfile 
if (defined($exptype) && ($exptype eq "continuous")) {
  $expfile = "$expfile\_PAGE/outquantized.txt";
}

$ta->loadFile($expfile);
$a_ref = $ta->getArray();

foreach my $r (@$a_ref) {
  if (defined($genes{$r->[0]})) {
    if (!defined($bin) || ((@a_bins > 0) && (Sets::in_array($r->[1], @a_bins)))) {
      print "$r->[0]\t$h_ref_desc->{$r->[0]}->[1]\t$h_ref_desc->{$r->[0]}->[2]";
      print "\t$r->[1]\n"
    }
  }
}
