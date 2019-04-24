use lib "$ENV{FIREDIR}/SCRIPTS";



#
# input: summary, expfile, profiles, 
# 
#


use Table;
use Sets;

use strict;

if (@ARGV == 0) {
  die "Usage: perl generate_motif_occurrence_summary.pl -expfile FILE -expfile_q FILE --summaryfile FILE -profiles FILE -pvmatrix FILE -columnfile FILE  -distmatrix FILE -oriematrix FILE -gomotif FILE -goindex FILE -profiles_orth FILE -orthologs FILE\n";
}

my $progdir     = "$ENV{FIREDIR}/PROGRAMS";


my $pvmatrix    = Sets::get_parameter(\@ARGV, "-pvmatrix");
my $columnfile  = Sets::get_parameter(\@ARGV, "-columnfile");
my $expfile     = Sets::get_parameter(\@ARGV, "-expfile");

my $distmatrix  = Sets::get_parameter(\@ARGV, "-distmatrix");
my $oriematrix    = Sets::get_parameter(\@ARGV, "-oriematrix");

my $quantized   = Sets::get_parameter(\@ARGV, "-quantized");

my $expfile_q   = undef;
if (Sets::exist_parameter(\@ARGV, "-expfile_q") == 1) {
  $expfile_q  =  Sets::get_parameter(\@ARGV, "-expfile_q");
}

my $summaryfile = Sets::get_parameter(\@ARGV, "-summaryfile");
my $profile     = Sets::get_parameter(\@ARGV, "-profiles");
my $rna         = undef;
if (Sets::exist_parameter(\@ARGV, "-rna") == 1) {
  $rna         = Sets::get_parameter(\@ARGV, "-rna");
}
my $gomotifs = undef;
if (Sets::exist_parameter(\@ARGV, "-gomotifs") == 1) {
  $gomotifs     = Sets::get_parameter(\@ARGV, "-gomotifs");
}

my $profiles_orth = undef;
if (Sets::exist_parameter(\@ARGV, "-profiles_orth") == 1) {
  $profiles_orth     = Sets::get_parameter(\@ARGV, "-profiles_orth");
}


my $mapping_orth = undef;
if (Sets::exist_parameter(\@ARGV, "-mapping_orth") == 1) {
  $mapping_orth     = Sets::get_parameter(\@ARGV, "-mapping_orth");
}

my $orthologs = undef;
if (Sets::exist_parameter(\@ARGV, "-orthologs") == 1) {
  $orthologs     = Sets::get_parameter(\@ARGV, "-orthologs");
}

my $goindex = undef;
if (Sets::exist_parameter(\@ARGV, "-goindex") == 1) {
  $goindex = Sets::get_parameter(\@ARGV, "-goindex");
}

my $rootdir     = ".";

my $ta = Table->new;

#
# read in column file
#
$ta->loadFile($columnfile);
my $a_ref_cols = $ta->getColumn(0);
my %COL = ();
for (my $i=0; $i<@$a_ref_cols; $i++) {
  $COL{ $a_ref_cols->[$i] } = $i;
}

#
# read in pvmatrix
#
$ta->loadFile($pvmatrix);
my $a_ref_pv = $ta->getArray();
shift @$a_ref_pv;
my %PV = ();
foreach my $r (@$a_ref_pv) {
  my $m = shift @$r;
  $PV{ $m } = $r;
}

#
# read in distmatrix
#
$ta->loadFile($distmatrix);
my $a_ref_distpv = $ta->getArray();

my %DISTPV = ();
foreach my $r (@$a_ref_distpv) {
  my $m = shift @$r;
  push @{ $DISTPV{ $m } }, $r;
}

#
# read in oriematrix
#
$ta->loadFile($oriematrix);
my $a_ref_oriepv = $ta->getArray();

my %ORIEPV = ();
foreach my $r (@$a_ref_oriepv) {
  my $m = shift @$r;
  push @{ $ORIEPV{ $m } }, $r;
}


#
# read in summary
#
$ta->loadFile($summaryfile);
my $a_ref_mo = $ta->getArray();
my @MOTIFS    = ();
my %CLUSTERS  = ();
foreach my $r (@$a_ref_mo) {
  push @MOTIFS, $r->[0];
  for (my $i=12; $i<@$r; $i++) {
    push @{ $CLUSTERS{ $r->[0] } }, $r->[$i];
  }	
}

#
# read in gomotifs
#
my $h_ref_go = undef;
if (defined($gomotifs)) {
  $ta->loadFile($gomotifs);
  $h_ref_go = $ta->getIndex(0);
}


# read in goindex
my $h_ref_goindex = undef;
if (defined($goindex)) {
  $ta->loadFile($goindex);
  $h_ref_goindex = $ta->getIndexShifted();
}


# read in expfile
$ta->loadFile($expfile);
my $h_ref_exp = $ta->getIndexKV(0,1);

# 
my $h_ref_exp_q = undef;
if (defined($expfile_q)) {
  $ta->loadFile($expfile_q);
  $h_ref_exp_q = $ta->getIndexKV(0,1);
}

# read in profiles
my %HITS = ();
$ta->loadFile($profile);
my $a_ref_pro = $ta->getArray();
foreach my $r (@$a_ref_pro) {

  if (defined($h_ref_exp->{$r->[1]})) {

    #print "defined $r->[1]\n";
  
    my $e = undef;
    if ($quantized == 1) {
      $e = $h_ref_exp->{$r->[1]}; 
    } else {
      $e = $h_ref_exp_q->{$r->[1]}; 
    }
    
    #print "e=$e\n";

    my $p = 0;
    #if (Sets::in_array($e, @{ $CLUSTERS{ $r->[0] } })) { 
    my $pv = $PV{ $r->[0] }[ $COL{$e} ];
    
    #
    # determine position bias pv
    #

    my $pvd = 0;
    foreach my $s (@{ $DISTPV{$r->[0] }}) {
      
      my ($i1, $i2) = $s->[1] =~ /\[([\d\.]+)\;([\d\.]+)\]/;
      die "Pb with parsing $s->[1]" if (!defined($i1));

      #if ($rna == 0) {
      if (($s->[0] == $e) && ($r->[2] >= $i1) && ($r->[2] <= $i2)) {
	$pvd = $s->[2];
      }
      #}
    }

    #
    # determine orientation bias pv
    #
    my $pvo = 0;
    foreach my $s (@{ $ORIEPV{$r->[0] }}) {
      #if ($rna == 0) {
	if ($s->[0] == $e) {
	  if ($r->[3] == 1) {
	    $pvo = $s->[1];
	  } else {
	    $pvo = $s->[2];
	  }
	}
      #}
    }
    
    my $comp = $pv + $pvd + $pvo;

    my @a_tmp = ($r->[1], $h_ref_exp->{$r->[1]}, $e, $comp, $pv, $pvd, $pvo, ($rna==0?-$r->[2]:$r->[2]), $r->[4]);
    push @{ $HITS{ $r->[0] }{ $r->[1] } }, \@a_tmp;
   
    #}
   
  } else {
    #print "$r->[1] (" . $h_ref_exp->{$r->[1]} . ") no defined\n";
  }
}


my %ORTH_HITS = ();
my %ORTH_PRES = ();

if (defined($profiles_orth)) {
  
  # 1. get gene correspondence
  my $h_ref_map = undef;
  if (defined($mapping_orth)) {
    $ta->loadFile($mapping_orth);
    $h_ref_map = $ta->getIndexKV(1,0);
  }
  

  # 2. read in profiles
  $ta->loadFile($profiles_orth);
  my $a_ref_orth = $ta->getArray();
  foreach my $r (@$a_ref_orth) {
    my $re = shift @$r;  
    
    my $g = $r->[0];

    if (defined($h_ref_map) && defined($h_ref_map->{$r->[0]})) {
      $g = $h_ref_map->{$r->[0]}; 
    }

    push @{ $ORTH_HITS{ $re }{ $g } }, $r;

    $ORTH_PRES{$g} = 1;
    
  }

}



# 
# 
#
my @HITS = ();
foreach my $re (@MOTIFS) {
  print "$re\n\n";
  
  print "Best GO enrichment: " . $h_ref_go->{$re}->[1] . "\n\n";
  
  print "Genes/ORF that have the motifs (genes are sorted by max(pa+pd+po)):\n\n";

  my @AR = ();
  foreach my $g (sort(keys(%{ $HITS{$re} }))) {
    push @AR, @{$HITS{$re}{$g}};
  }
  
  @AR = sort { $b->[3] <=> $a->[3] } @AR;  # sort sites based on composite score

  my @GENES = ();
  foreach my $g (@AR) {
    push @GENES, $g->[0] if (!Sets::in_array($g->[0], @GENES));
  }

  foreach my $o (@GENES) {
    #if ($quantized == 1) {

    print "$o\n";
    if (defined($h_ref_goindex) && defined($h_ref_goindex->{$o}) && Sets::in_array($h_ref_go->{$re}->[2], @{$h_ref_goindex->{$o}})) {
      my $txt  = $h_ref_go->{$re}->[1];
      $txt =~ s/\,\ p.+$//;
      print "  Has best GO annotation: " . $txt . "\n";
    }
    print "  Expression: " . $HITS{ $re }{ $o }->[0]->[1] . "\n";
    print "\n";
    print sprintf("%10s\t%3s\t%3s\t%3s\n", "Position", "pa", "pd", "po");
    foreach my $g (@{ $HITS{ $re }{ $o } }) {
      print sprintf("%10d\t%2.1f\t%2.1f\t%2.1f\t%s", $g->[7], $g->[4], $g->[5], $g->[6], $g->[8]);
      print "\n";
    }
    print "\n";

    if (defined($profiles_orth)) {

      print "  Ortholog: ";
      if (defined($ORTH_PRES{$o})) {
	print "$o";
      } else {
	print "none";
      }
      
      print "\n";

      my $popo =  $ORTH_HITS{ $re }{ $o };
      
      foreach my $po (@$popo) {
	print sprintf("%10d\t   \t   \t   \t%s\n", ($rna==0?-$po->[1]:$po->[1]), $po->[3]);
      }

      print "\n";
    }

  }
  print "\n";
  print "---------------------------------------------------------------------------\n\n";
}
