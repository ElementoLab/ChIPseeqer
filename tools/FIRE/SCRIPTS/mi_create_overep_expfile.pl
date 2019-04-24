use lib "$ENV{FIREDIR}/SCRIPTS";

use Table;
use Sets;

use strict;

if (@ARGV == 0) {
  die "Usage: perl mi_create_overep_expfile.pl -repfile FILE -expfile FILE -outexpfile FILE -quantized INT \n";
}



my $progdir = "$ENV{FIREDIR}/PROGRAMS";

my $repfile   = undef;
if (Sets::exist_parameter(\@ARGV, "-repfile") == 1) {
  $repfile    = Sets::get_parameter(\@ARGV, "-repfile");
}

my $expfile   = undef;
if (Sets::exist_parameter(\@ARGV, "-expfile") == 1) {
  $expfile   = Sets::get_parameter(\@ARGV, "-expfile");
}

my $outexpfile     = undef;
if (Sets::exist_parameter(\@ARGV, "-outexpfile") == 1) {
  $outexpfile = Sets::get_parameter(\@ARGV, "-outexpfile");
}

my $quantized      = Sets::get_parameter(\@ARGV, "-quantized");




#  parse MOTIF REP files !
#

my %pvalues   = ();
my %densities = ();
my @labels    = ();
my %ebins     = ();
my @MOTIFS    = ();

parse_motif_rep_file($repfile, $quantized, \%pvalues, \%densities, \@labels, \%ebins, \@MOTIFS);

my $ta = Table->new;
$ta->loadFile($expfile);
my $a_ref = $ta->getArray();
my $r = shift @$a_ref;

open TUO, ">$outexpfile" or die "Cannot open $outexpfile.\n";

print TUO "ORF";
foreach my $m (@MOTIFS) {
  print TUO "\t$m";
}
print TUO "\n";
foreach my $r (@$a_ref) {
  print TUO "$r->[0]";
  foreach my $m (@MOTIFS) {
    if (Sets::in_array($r->[1], @{ $ebins{$m} })) {
      print TUO "\t1";
    } else {
      print TUO "\t0";
    }
  }
  print TUO "\n";
}



sub parse_motif_rep_file {
  
  my ($repfile, $quantized, $h_ref_pvalues, $h_ref_densities, $a_ref_labels, $h_ref_ebins, $a_ref_motifs) = @_;



  my @MOTIFS      = ();
  my $cnt         = 1;
  my %FREQ        = ();
  my @N           = ();
  my %avg_density = ();
  my %sample_size = ();

  #
  # read in motif repfile
  #
  
  my $ta = Table->new;
  
  $ta->loadFile($repfile);
  my $a_ref = $ta->getArray();
  
  foreach my $r (@$a_ref) {
    
    #
    #  motif is always first column
    #
    push @MOTIFS, $r->[0] if !Sets::in_array($r->[0], @MOTIFS);
    
    if ($quantized == 1) {
      push @{ $FREQ{ $r->[0] } }, $r->[4];
      push @N, $r->[5];  # number of item per cluster
      
      $sample_size{ $r->[0] } += $r->[5];
      $avg_density{ $r->[0] } += int( 0.5 + $r->[4] * $r->[5] );
      
    } else {
      push @{ $FREQ{ $r->[0] } }, $r->[6];
      push @N, $r->[7];  # number of item per cluster
      
      $sample_size{ $r->[0] } += $r->[7];
      $avg_density{ $r->[0] } += int( 0.5 + $r->[6] * $r->[7] );
    }
    
    my $label = undef;
    if ($quantized == 1) {
      $label = "C$r->[1]";
    } else {
      $label = "[$r->[4];$r->[5]]";
    }
    
    push @$a_ref_labels, $label if (!Sets::in_array($label, @$a_ref_labels));

    $cnt ++;
  }
  

  #
  # calculate average density
  #
  foreach my $k (keys(%avg_density)) {
    $avg_density{$k} = $avg_density{$k} / $sample_size{$k};
  }
  
  
  
  #
  # record significant exp bins for each motif
  #

  
  my $cnt = 0;
  foreach my $k (@MOTIFS) {
    
    #print "Processing $k.                        \n";
    
    my $n = scalar( @{ $FREQ{$k} } ); 
    
    #print "I have $n bins\n";

    for (my $i=0; $i<$n; $i++) {
      
      

      my $bk = undef;
      my $bn = undef;
      
      my $the_ebin = undef;
      
      $bk = int( 0.5 + $FREQ{$k}->[$i] * $N[$i] );
      $bn = $N[$i];
      $the_ebin = $i;
      
      my $bN = $avg_density{$k};    
      my $p2 = undef;
      
      my $p1 = binom_test_greater( $bk, $bn, $bN, \$p2);
      
      my $p  = undef;
      
      if ($p1 < 1e-100) {
	$p1 = 1e-100;
      }
      if ($p2 < 1e-100) {
	$p2 = 1e-100;
      }
      if ($p1 < $p2) {
	$p = sprintf("%4.3f", - Sets::log10($p1));
      } else {
	$p = sprintf("%4.3f", + Sets::log10($p2));
      }
      
      if ($p1 < 0.05 / $n) {
	push @{ $h_ref_ebins->{ $k } }, $the_ebin;
      }
      
      $h_ref_pvalues   ->{ $k }->[ $i ] = $p;

      $h_ref_densities ->{ $k }->[ $i ] = $FREQ{$k}->[$i];
      
    }
    
    
    $cnt ++;    
 
    #print "Done\n";
  }

  @$a_ref_motifs = @MOTIFS;
  
}



sub binom_test_greater {
  
  my ($n, $k, $p, $p1) = @_;

  my $todo = "$progdir/binom_test_greater $n $k $p 1";
  my $out = `$todo`; 
  $out =~ s/[\n\r]//g;

  my @a = split /\t/, $out, -1;
  
  $out = $a[0];

  $$p1 = $a[1];
  
  return $out;

}



sub parse_distance_file {
  my ($distfile, $RDIST, $RORIE, $RCOPY) = @_;
  
  my $o_t = 1;

  #
  #  load distance info
  #
  $ta->loadFile($distfile);
  my $a_ref_dist = $ta->getArray();

  my %TMPDIST   = ();
  my @MOTIFS    = ();
  foreach my $r (@$a_ref_dist) {
    my $mo = shift @$r; 
    push @MOTIFS, $mo if (!Sets::in_array($mo, @MOTIFS));
    my $ke = shift @$r;
    $TMPDIST{$mo}{$ke} = $r;
  }
  
  foreach my $m (@MOTIFS) {
    
    #if (($TMPDIST{ $m }{"d_med"}->[2] <= 10) || ($TMPDIST{ $m }{"d_avg"}->[2] <= 10)) {

    if ($TMPDIST{ $m }{"d_avg"}->[2] <= 100) {    
      $RDIST->{ $m } = 1;
    } else {
      $RDIST->{ $m } = 0;
    }
    
    if (($TMPDIST{ $m }{"o5"}->[2] == 0 ) && ($TMPDIST{ $m }{"o3"}->[2] >= $o_t)) {
      $RORIE->{ $m } = 1;
    } elsif (($TMPDIST{ $m }{"o5"}->[2] >= $o_t) && ($TMPDIST{ $m }{"o3"}->[2] ==  0)) {
      $RORIE->{ $m } = 2;
    } else {
      $RORIE->{ $m } = 0;
    }
    
    my $i=1; my $bi = undef; my $zz = -10000; my $mm = -10000;
    while (defined($TMPDIST{ $m }{"c$i"})) {
      if (($TMPDIST{ $m }{"c$i"}->[3] > $zz) && ($TMPDIST{ $m }{"c$i"}->[1] > $mm))  {
	$bi = $i; $mm = $TMPDIST{ $m }{"c$i"}->[1]; $zz = $TMPDIST{ $m }{"c$i"}->[3];
      }
      $i++;
    }
    
    $RCOPY->{$m} = $bi;
  }
  
}
