use lib "$ENV{FIREDIR}/SCRIPTS";


use Table;
use Sets;
use Getopt::Long;
use strict;

my $namefile    = undef;
my $acefiles    = undef;
my $summaryfile = undef;

my $progdir     = "$ENV{FIREDIR}/PROGRAMS";

my $compareace  = "$progdir/MyCompareAce";
die "Please make sure that MyCompareAce is compiled.\n" if ((! -e $compareace) && (! -e "$compareace.exe"));

my $rna         = undef;
my $shufflecols = 0;
my $outfile     = undef;
my $threshold   = 0.8;
my $shufflemotifs = undef;

#srand(12345);

if (@ARGV == 0) {
  die "perl annotate_motifs_using_ace_library.pl --acefiles --summaryfile --namefile= --shufflecols= --rna=0 \n";
}

GetOptions ('acefiles=s'        => \$acefiles,
	    'summaryfile=s'     => \$summaryfile,
	    'namefile=s'        => \$namefile,
	    'shufflecols=s'     => \$shufflecols,
	    'shufflemotifs=s'   => \$shufflemotifs,
	    'outfile=s'         => \$outfile,
	    'threshold=s'       => \$threshold,
	    'rna=s'             => \$rna);



die "Please define -rna\n" if (!defined($rna));

my $a_ref_library = Sets::getFiles($acefiles);


my $ta = Table->new;



die "Please enter summaryfile\n" if (! -e $summaryfile);
$ta->loadFile($summaryfile);
my $a_ref_sum = $ta->getArray();
my $a_ref_mot = $ta->getColumn(0); 
my $cnt = 1;

if (defined($shufflemotifs)) {
  
  my @tmpmot = ();
  foreach my $r (@$a_ref_sum) {
    my $re = $r->[0];	   
    my $more = $re;
    $more =~ s/^\.+//;
    $more =~ s/\.+$//;
    
    my $a_ref_re     = Sets::get_array_from_re($more);
    my $a_ref_re_shu = Sets::shuffle_array($a_ref_re);
    
    $r->[0] = join("", @$a_ref_re_shu);
    
    push @tmpmot, $r->[0];
    
  }

  @$a_ref_mot = @tmpmot;
}


my $h_ref_names = {};
if (defined($namefile)) {
  $ta->loadFile($namefile);
  $h_ref_names = $ta->getIndexKV(0,1);
} else {
  foreach my $r (@$a_ref_mot) {
    $h_ref_names->{ $r } = "-";
  }
}

#foreach my $r (keys(%$h_ref_names)) {
#  print "$r\n";
#}

foreach my $r (@$a_ref_sum) {
    

  next if (($rna != 2) && ($r->[1] != $rna));  

  my $re = $r->[0];	   
  my $more = $re;
  $more =~ s/^\.+//;
  $more =~ s/\.+$//;

  # output WMed motif to tmp file
  my $mo = Sets::myre2scanace($more);
  my $mn = Sets::getTempFile("/tmp/mytmpmotif");
  open OUT, ">$mn" or die "cannot open mot file for output.\n";
  print OUT "Motif 1\n$mo";
  close OUT;
  
  my $best_c = -10000;
  my $best_m = undef;
  foreach my $m (@$a_ref_library) {
    
    #print "$m\n";
    
    my $mym = $m;
    if ($shufflecols == 1) {
      my $txt = readAndShuffleAceMotif($m);
      open OUT, ">SHU";
      print OUT "Motif 1\n$txt";
      close OUT;
      $mym = "SHU";
    }

    my $s_todo = "$compareace $mn \"$mym\" -simple ";
    if ($rna == 1) {
      $s_todo .= " -ss";
    }
    my $score = `$s_todo`;
    chomp $score;
    if ($score > $best_c) {
      $best_c = $score;
      $best_m = $m;
    }
  }
    
  unlink $mn;

  my $momo = Sets::filename($best_m);
  $momo =~ s/\.ace//g;
  
  if ($best_c >= $threshold) {
    
    
    if ($h_ref_names->{ $re } eq "-") {
      $h_ref_names->{ $re }  = $momo;
    } else {
      $h_ref_names->{ $re } .= $momo;
    }	
  } else {
    #print "$re\t-\n";
  }
  

  $cnt ++;

  
}



open OO, ">$outfile" or die "Cannot open $outfile.\n";
foreach my $r (keys(%$h_ref_names)) {
  print OO "$r\t$h_ref_names->{$r}\n";
  print    "$r\t$h_ref_names->{$r}\n";
}
close OO;



sub readAndShuffleAceMotif {
  my ($f) = @_;

  my $txt = "";

  open IN, $f or die "Cannot open $f.\n";
  my @a = <IN>; chomp @a;
  shift @a;

  my $n = scalar(@a);
  my $m = undef;

  my @H = ();
  foreach my $s (@a) {
    my @b = split //, $s;
    $m = scalar(@b); #print "m=$m\n";
    for (my $i=0; $i<scalar(@b); $i++) {
      
      push @{ $H[$i] }, $b[$i];
    }
  }
  close IN;

  my $a_ref = Sets::shuffle_array(\@H);
  
  for (my $i=0; $i<$n; $i++) {
    
    for (my $j=0; $j<$m; $j++) {
      $txt .= $a_ref->[$j]->[$i];
    }
    $txt .=  "\n";
    
  }
  
  return $txt;
}
