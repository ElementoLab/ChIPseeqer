#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Getopt::Long;
use Table;
use strict;

my $todo  = "";

my $dodef      = 1;
if (grep(/--donothing=1/, @ARGV)) {
  $dodef = 0;
  @ARGV = grep !/--donothing/, @ARGV;
}

if (@ARGV == 0) {
  die "Usage: meta-page.pl [same args as PAGE]\n";
}

my $expfile    = undef;
my $dogo       = $dodef;
my $doint      = $dodef;
my $dorep      = $dodef;
my $dokegg     = $dodef;
my $dobiocarta = $dodef;
my $docmap     = $dodef;
my $docombine  = 0;
my $dosystem   = 1;
my $ebins      = undef;
my $verbose    = 0;
my $exptype    = 1;
my $randomize  = undef;
my $max_p      = undef;

my @myARGV = @ARGV;

GetOptions("expfile=s"    => \$expfile,
	   "exptype=s"    => \$exptype,
	   "dokegg=s"     => \$dokegg,
	   "docombine=s"  => \$docombine,
	   "dobiocarta=s" => \$dobiocarta,
	   "max_p=s"      => \$max_p,
	   "dogo=s"       => \$dogo,
	   "verbose=s"    => \$verbose,
	   "ebins=s"      => \$ebins,
	   "dosystem=s"   => \$dosystem,
	   "doint=s"      => \$doint,
	   "dorep=s"      => \$dorep,
	   "docmap=s"     => \$docmap,
	   "randomize=s"  => \$randomize);

@myARGV = grep !/docmap/, @myARGV;
push @myARGV, "--suffix=pathways";
push @myARGV, "--outprofiles=1";

my $args = join("\t", @myARGV);

my @pathway_list = ();

# GO
if ($dogo == 1) {

  print "GO Biological processes.\n";

  # P
  $todo = "page.pl $args --pathways=human_go_orf";
  push @pathway_list, "human_go_orf";
  mysystem($todo);

  # C
  print "GO Cellular components.\n";
  $todo = "page.pl $args --pathways=human_go_orf --cattypes=C";
  push @pathway_list, "human_go_orf";
  mysystem($todo);
  
  # F
  print "GO Biochemical functions.\n";
  $todo = "page.pl $args --pathways=human_go_orf --cattypes=F";
  push @pathway_list, "human_go_orf";
  mysystem($todo);

}


# interactions
if ($doint == 1) {
  
  print "HPRD interactions.\n";
  $todo = "page.pl $args --pathways=HPRD_interactions";
  push @pathway_list, "HPRD_interactions";
  mysystem($todo);
}

# repeats
if ($dorep == 1) {
  
  print "Repeats (ORF level)\n";
  $todo = "page.pl $args --pathways=human_repeats_orf";
  push @pathway_list, "human_repeats_orf";
  mysystem($todo);
}


# KEGG
if ($dokegg == 1) {
  print "KEGG pathways.\n";
  $todo = "page.pl $args --pathways=kegg";
  push @pathway_list, "kegg";
  mysystem($todo);
}


# Biocarta
if ($dobiocarta == 1) {
  print "Biocarta pathways.\n";
  $todo = "page.pl $args --pathways=biocarta";
  push @pathway_list, "biocarta";
  mysystem($todo);
}

# connnectivity map

if ($docmap == 1) {
  print "Connectivity map.\n";
  $todo = "e-page.pl $args --profiles=cmap";
  push @pathway_list, "cmap";
  mysystem($todo);
}


if ($docombine == 1) {

  # build master-table
  
  
  my $ta = Table->new;
  my @matrices = ();
  my @HEADER   = ();
  my @HEADER_n = ();
  my %genes    = ();
  
  foreach my $p (@pathway_list) {
    my $f = "$expfile.$p\_PAGE/signif.profiles.txt";
    
    $ta->loadFile($f);
    my $a_ref = $ta->getArray();
    
    # add header
    my $h = shift @$a_ref;
    shift @$h;
    
    foreach my $c (@$h) {
      if ($p eq "cmap") {
	$c = "C $c";
      } else {
	$c = "D $c";
      }
  }
    
    push @HEADER, @$h;
    push @HEADER_n, scalar(@$h);
    
    my $i = 0;
    my $h_matrix = {};
    foreach my $r (@$a_ref) {
      my $n = shift @$r;
      $genes{ $n } = 1;
      $h_matrix->{$n} = $r;
      $i ++;
    }
    
    push @matrices, $h_matrix;
  }
  
  
  my %MATRIX = ();
  
  foreach my $g (keys(%genes)) {
    
    # traverse matrices
    my $k = 0;
    foreach my $h (@matrices) {
      if (defined($h->{$g})) {
	push @{ $MATRIX{ $g } }, @{ $h->{$g} };
      } else {
	for (my $i=0; $i<$HEADER_n[$k]; $i++) {
	  push @{ $MATRIX{ $g } }, "NA";
	}	
      }
      $k++;
    }
    
  }
  
  
  
  mkdir "$expfile\_META" if (! -e "$expfile\_META");
  open OUT, ">$expfile\_META/combined.matrix.txt";
  
  print OUT "GENE\t" . join("\t", @HEADER) . "\n";
  foreach my $g (keys(%genes)) {
    print OUT "$g\t" . join("\t", @{ $MATRIX{$g} }) . "\n";
  }
  close OUT;
  
  # calculate all pairwise interactions
  $todo = "$ENV{PAGEDIR}/PROGRAMS/mi_analyze_matrix $expfile\_META/combined.matrix.txt > $expfile\_META/combined.matrix.interactions.txt";
  system($todo);
  
  # draw
  $todo = "$ENV{PAGEDIR}/SCRIPTS/draw_profile_interaction_matrix.pl $expfile\_META/combined.matrix.txt > $expfile\_META/combined.matrix.interactions.txt";
  system($todo);
  
  
  sub mysystem {
    my ($c) = @_;
    
    if ($dosystem == 1) {
      return system($c);
    } else {
      return -1;
    }	
  }
  
}
