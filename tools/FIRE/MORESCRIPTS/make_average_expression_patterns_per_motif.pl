use lib "$ENV{FIREDIR}/SCRIPTS";

#
# 1. get output of micombine as input
# 2. do a GO analysis on them
# 3
#

use Sets;
use GroupEnrichment;
use Getopt::Long;
use Table;

use strict;


my $re             = undef;
my $fastafile      = undef;
my $goindex        = undef;
my $gonames        = undef;
my $maxcatsize     = 500;
my $pmin           = 0.01;
my $outgo          = undef;
my $outfile        = undef;

my $profiles       = undef;
my $expfile        = undef;

my $clusterfile    = undef;
my $usemodule      = undef;
my $datafile       = undef;
my $shuffle        = 0;
my $where          = 'DNA_RNA';

if (@ARGV == 0) {
  die "Usage: perl make_average_expression_patterns_per_motif.pl --summaryfile=FILE --expfile=FILE --profiles=FILE --goindex=FILE --gonames=FILE\n";
}

GetOptions (
	    'expfile=s'        => \$expfile,
	    'datafile=s'       => \$datafile,
	    'shuffle=s'        => \$shuffle,
	    'where=s'          => \$where
	    );


use Table;

my $ta = Table->new;

my $expfile_f    = Sets::filename($expfile);
my $summaryfile  = "$expfile\_FIRE/$where/$expfile.summary";
my $profiles_dna = "$expfile\_FIRE/DNA/$expfile.profiles";
my $profiles_rna = "$expfile\_FIRE/RNA/$expfile.profiles";

my $targetdir    = "$expfile\_FIRE/$where/TARGETS";
mkdir $targetdir if (! -e $targetdir);

#
# read in clusters associated to each motif
#


my %CLU = ();
my @MOTIFS = ();
$ta->loadFile($summaryfile);
my $a_ref_sum = $ta->getArray();
my $hasdna = 0;
my $hasrna = 0;
foreach my $r (@$a_ref_sum) {
  push @MOTIFS, $r->[0];
  $hasrna = 1 if ($r->[0] == 1);
  $hasdna = 1 if ($r->[0] == 0);
  for (my $i=12; $i<@$r; $i++) {
    $CLU{$r->[0]}{$r->[$i]} = 1;
  }
}

#
#  read in expfile
#
$ta->loadFile($expfile);
my $h_ref_exp = undef;

if ($shuffle == 0) {
  $h_ref_exp = $ta->getIndexKV(0,1);
} else {

  my $a_ref_col0 = $ta->getColumn(0);
  my $a_ref_col1 = $ta->getColumn(1);
  my $h1         = shift @$a_ref_col0;
  my $h2         = shift @$a_ref_col1;

  my $a_ref_shu1 = Sets::shuffle_array($a_ref_col1);

  for (my $i=0; $i<@$a_ref_col0; $i++) {
    die "dying ..\n" if (!defined($a_ref_shu1->[$i]));
    $h_ref_exp->{ $a_ref_col0->[$i] } = $a_ref_shu1->[$i];
  }
}

#my $a_ref_exp = $ta->getColumn(0);

my %CLU_GENES = ();



if ($hasdna == 1) {

  $ta->loadFile($profiles_dna);
  my $a_ref_pro = $ta->getArray();
  foreach my $r (@$a_ref_pro) {
    # get the expression cluster
    my $c = $h_ref_exp->{ $r->[1] };
    next if (!defined($c));
    
    # motif
    my $m = $r->[0];
    
    # if this cluster that this gene belongs to is a cluster where the motif is enriched, keep gene 
    if (defined($CLU{$m}{$c})) {
      push @{ $CLU_GENES{ "$m" }}, $r->[1];
    }
  }

}

#
# load RNA profiles
#

if ($hasrna == 1) {
  $ta->loadFile($profiles_rna);
  my $a_ref_pro = $ta->getArray();
  foreach my $r (@$a_ref_pro) {
    # get the expression cluster
    my $c = $h_ref_exp->{ $r->[1] };
    next if (!defined($c));
    
    # motif
  my $m = $r->[0];
    
    # if this cluster that this gene belongs to is a cluster where the motif is enriched, keep gene 
    if (defined($CLU{$m}{$c})) {
      push @{ $CLU_GENES{ "$m-" }}, $r->[1];
    }
  }
  
}

#
# traverse motifs, remove duplicates, calculate average expression patterns
#

#my $avgfile = "$targetdir/avgexpfile";
#open OUT, ">$avgfile" or die "Cannot open $avgfile";

my $h_ref_names = {};
if (-e "$datafile.names") {
  my $cnt = 0;
  open IN, "$datafile.names" or die "dkjfkje\n";
  while (my $l = <IN>) {
    chomp $l;
    if (!defined($h_ref_names->{$l})) {
      $h_ref_names->{$l} = $cnt;
    }
    $cnt++;
  }
  close IN;
}



foreach my $k (keys(%CLU_GENES)) {

  $CLU_GENES{$k} = Sets::removeDuplicates($CLU_GENES{$k});

  print "Motif $k, found " . scalar( @{ $CLU_GENES{$k} } ) . " unique targets.\n"; 

  if (scalar( @{ $CLU_GENES{$k} } ) > 0 ) {

    print "Saving ..\n";

    my $ff = "$expfile\_FIRE/$where/TARGETS/$k.txt";
    Sets::writeSet($CLU_GENES{$k}, $ff);
    my $h_ref_H = Sets::getIndexFromArrayRef($CLU_GENES{$k});

    #
    # calculate average expression pattern (reopen datafile each time, this is a long process)
    #

    my @avg = ();
    my $cnt = 0;

    if (! -e "$datafile.idx") {

      open IN, $datafile or die "Cannot open $datafile.\n";

      while (my $l = <IN>) {

	my ($n) = $l =~ /^(.+?)\t/;
	next if ($n eq 'GENE');

	if (defined($h_ref_H->{$n})) {

	  chomp $l;
	  my @a = split /\t/, $l, -1;
	

	  shift @a;
	  if ($cnt == 0) {
	    @avg = @a;
	  } else {
	    Sets::addArray2ToArray1(\@avg, \@a);
	  }
	  $cnt ++;
	}	
      }
      close IN;

    } else {

      open(FILE, "< $datafile")      or die "Can't open $datafile for reading: $!\n";
      open(INDEX, "< $datafile.idx") or die "Can't open $datafile.idx for reading: $!\n";

      print "Using $datafile.idx.\n";
      
      foreach my $n (@{ $CLU_GENES{$k} }) {
	my $l = line_with_index(*FILE, *INDEX, $h_ref_names->{$n});
	chomp $l;
	
	my @a = split /\t/, $l, -1;

	print "Got line for $n ($h_ref_names->{$n}) .. ";
	print "Line has " . scalar(@a) . " elements, starts with $a[0], $a[1]...\n";


	shift @a;
	if ($cnt == 0) {
	  @avg = @a;
	} else {
	  Sets::addArray2ToArray1(\@avg, \@a);
	}
	$cnt ++;
      }

    }

    print "Found $cnt targets.\n";

    my $ft = "$k.txt.avgexp";
    $ft =~ s/^\.+//;
    
    $ft = "$expfile\_FIRE/$where/TARGETS/$ft";
    if ($shuffle == 1) {
      $ft .= ".shu";
    }
    open OUT, ">$ft" or die "Cannot open $ft.\n";
    
    for (my $i=0; $i<@avg; $i++) {
      $avg[$i] = sprintf("%4.3f", $avg[$i] / $cnt);
      print OUT $avg[$i] . "\n";
    }
    close OUT;

    printf "Wrote $ft.\n";
  
  } # end (if size targets > 0)

}



sub line_with_index {
  my $data_file   = shift;
  my $index_file  = shift;
  my $line_number = shift;
  
  my $size;               # size of an index entry
  my $i_offset;           # offset into the index of the entry
  my $entry;              # index entry
  my $d_offset;           # offset into the data file
  
  $size = length(pack("L!", 0));
  $i_offset = $size * $line_number;
  seek($index_file, $i_offset, 0) or return;
  read($index_file, $entry, $size);
  $d_offset = unpack("L!", $entry);
  seek($data_file, $d_offset, 0);
  return scalar(<$data_file>);
}



