# usage 
#
# f=FIRE_DATA/HUMAN/EXPFILES/TISSUE/human_tissue_IclustPos.txt_FIRE/DNA_RNA/human_tissue_IclustPos.txt
#
#
# perl  mi_restrict_files_for_figures.pl  --matfile=$f.matrix --densityfile=$f.densities --clustfile=$f.clusters --columnsfile=$f.columns --restrict=restrict.human
# 
# perl fire.pl --expfiles=FIRE_DATA/HUMAN/EXPFILES/TISSUE/human_tissue_IclustPos.txt --exptype=discrete --species=human --dodef=0 --domicombine=2 --domidrawmatrix=1 --dodna=0 --dorna=0
# 
# perl  mi_restrict_files_for_figures.pl  --matfile=$f.matrix --densityfile=$f.densities --clustfile=$f.clusters --columnsfile=$f.columns --restrict=restrict.human --revert=1
#
#


# p-value
# densities
# columns file
# clusters

# retrict


use lib qw(PostScript-Simple-0.07/lib);
use Table;
use Sets;
use Getopt::Long;
use PostScript::Simple;
use strict;



my $matfile        = undef;
my $densityfile    = undef;
my $clustfile      = undef;
my $columnsfile    = undef;
my $restrict       = undef;
my $revert         = undef;
my $samecolumns    = 0;

if (@ARGV == 0) {
  die "Usage: perl mi_restrict_files_for_figures.pl  --matfile=FILE --densityfile=FILE --clustfile=FILE --columnsfile=FILE --restrict=FILE --revert=INT --samecolumns=INT\n";
}

GetOptions ('matfile=s'        => \$matfile,
	    'densityfile=s'    => \$densityfile,
	    'restrict=s'       => \$restrict,
	    'clustfile=s'      => \$clustfile,
	    'columnsfile=s'    => \$columnsfile,
	    'samecolumns=s'    => \$samecolumns,
	    'revert=s'         => \$revert);


if ($revert == 1) {
  
  print "reverting to current .befrestrict files ...";
  #
  #  backup the 4 files if not yet dome
  #
  if ( -e "$matfile.befrestrict") {
    system("cp $matfile.befrestrict $matfile");
  }
  if ( -e "$densityfile.befrestrict") {
    system("cp $densityfile.befrestrict $densityfile");
  }
  if ( -e "$clustfile.befrestrict") {
    system("cp $clustfile.befrestrict $clustfile");
  }
  if ( -e "$columnsfile.befrestrict") {
    system("cp $columnsfile.befrestrict $columnsfile");
  }

  print "Done (hopefully).\n";
  
  exit;

}

#
#  backup the 4 files if not yet dome
#
if (! -e "$matfile.befrestrict") {
  system("cp $matfile $matfile.befrestrict");
}
if (! -e "$densityfile.befrestrict") {
  system("cp $densityfile $densityfile.befrestrict");
}
if (! -e "$clustfile.befrestrict") {
  system("cp $clustfile $clustfile.befrestrict");
}
if (! -e "$columnsfile.befrestrict") {
  system("cp $columnsfile $columnsfile.befrestrict");
}


my $ta = Table->new;

#
#  read in the restriction file
#
my $h_ref_rest = Sets::getIndex($restrict);

#
#  read in the columns file
#
my $a_ref_cols = Sets::readSet($columnsfile);

#
#  read in clusterfile
#
$ta->loadFile($clustfile);
my $a_ref_clust = $ta->getArray();

#
#  read in the matrix file
#
$ta->loadFile($matfile);
my $a_ref_M      = $ta->getArray();
my @newmat       = ();
my $h            = shift @$a_ref_M;
push @newmat, $h;

my @cols = (); 
foreach my $s (@cols) { 
  push @cols, 0; 
}

#
# find the minimum abs() in matrix
#
my $minabs = 100000;

foreach my $r (@$a_ref_M) {
  my $maxinrow = 0;
  for (my $i=1; $i<@$r; $i++) { 
    $maxinrow = abs($r->[$i]) if (abs($r->[$i]) > $maxinrow); 
  }
  if ($maxinrow < $minabs) {
    $minabs = $maxinrow;
  }
}

foreach my $r (@$a_ref_M) {
  if (defined($h_ref_rest->{$r->[0]})) {
    push @newmat, $r;
    for (my $i=1; $i<@$r; $i++) { 
      $cols[$i] = 1 if (abs($r->[$i]) > $minabs); 
    }
  }
}


#
# use @cols to keep only the good columns
#

if ($samecolumns == 0) {
  for (my $j=0; $j<@newmat; $j++) {
    my $r = $newmat[$j];
    my @a = (); push @a, $r->[0];
    for (my $i=1; $i<@$r; $i++) { 
      if ($cols[$i] == 1) {
	push @a, $r->[$i]; 
      }
    }
    $newmat[$j] = \@a;
  }
}

#
# ok, newmat contains the good stuff, now do the same for density
#
my @newdens = ();
$ta->loadFile($densityfile);
my $a_ref_D      = $ta->getArray();
my $h            = shift @$a_ref_D;
push @newdens, $h;

foreach my $r (@$a_ref_D) {
  if (defined($h_ref_rest->{$r->[0]})) {
    push @newdens, $r;
  }
}

if ($samecolumns == 0) {
  # keep the good cols
  for (my $j=0; $j<@newdens; $j++) {
    my $r = $newdens[$j];
    my @a = (); push @a, $r->[0];
    for (my $i=1; $i<@$r; $i++) { 
      if ($cols[$i] == 1) {
	push @a, $r->[$i]; 
      }
    }
    $newdens[$j] = \@a;
  }
}

#
# remove non-signif cols from columnsfile
#
my @newcols = @$a_ref_cols;

if ($samecolumns == 0) {
  @newcols = ();
  for (my $i=0; $i<@$a_ref_cols; $i++) {
    #print "$a_ref_cols->[$i] $cols[$i]\n";
    if ($cols[$i+1] == 1) {
      push @newcols, $a_ref_cols->[$i];
    }
  }
}

#
# remove motifs from clustfile
#
my @newclust = ();
foreach my $r (@$a_ref_clust) {
  if (defined($h_ref_rest->{$r->[0]})) {
    push @newclust, $r; 
  }
}


#
# and now, save to disk
#


#1 matrix
Sets::writeTable(\@newmat, "$matfile");

#2 density
Sets::writeTable(\@newdens, "$densityfile");

#3 columns
Sets::writeSet(\@newcols, "$columnsfile");

#4 clusters
Sets::writeTable(\@newclust, "$clustfile");



