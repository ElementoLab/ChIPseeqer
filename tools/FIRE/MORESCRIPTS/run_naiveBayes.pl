use lib "$ENV{FIREDIR}/SCRIPTS";
use Getopt::Long;
use Sets ;
use Data::Dumper ;

use Table ;
use strict;
use Subs ;
use Bayes ;

my $expfile = undef ;
my $profile = undef ;

GetOptions ('expfile=s'    => \$expfile,
	    'profile=s'    => \$profile) ;

my $namesfile = $profile ;
$namesfile =~ s/profiles$/motifnames/ ;
my $summaryfile = $profile ;
$summaryfile =~ s/profiles$/summary/ ;
my $distfile = $profile ;
$distfile =~ s/profiles$/dist.rep/ ;

my %motifnames ;
my %M_id ;

open I, "< $namesfile" or die "Couldn't: $namesfile";
print "Reading motif names...";
while(<I>){
    chomp;
    my ($motif, $name) = split(/\t/, $_);
    $motifnames{$motif} = $name ;
    my $cnt = scalar(keys %M_id) ;
    $M_id{$motif} = $cnt ;
}
close I;
print "Done\n";

open I, "< $profile" or die "Couldn't: $profile";
print "Reading motif profiles...";
my %status ;
while(<I>){
    chomp;
    my ($m, $g, $p, $o, @a) = split(/\t/, $_);
    $status{$m}{$g}{m} = 1 ;
    push(@{$status{$m}{$g}{p}}, $p) ;
    push(@{$status{$m}{$g}{o}}, $o) ;
}
close I;
print "Done\n";

open I, "< $summaryfile" or die "Couldn't: $summaryfile";
my @sig_clusts ;
my %pos_bias ;
my %ori_bias ;
while(<I>){
    chomp;
    my @a = split(/\t/, $_) ;
    $pos_bias{$a[0]} = $a[9] ;
    $ori_bias{$a[0]} = $a[10] ;
    for (my $i=12 ; $i<@a ; $i++){
	if (!(grep {$_ eq $a[$i]} (@sig_clusts))){
	    push(@sig_clusts, $a[$i]) ;
	}
    }
}

open I, "< $distfile" or die "Couldn't: $distfile";
while(<I>){
    chomp;
    my ($m, $l, $e, @a) = split(/\t/, $_) ;
    if ($l eq "d_avg" && $e eq "nan"){
	pop @a ;
	$status{$m}{bins} = \@a ;
    }
}

open I, "< $expfile" or die "Couldn't: $expfile" ;
print "Reading $expfile...";

open O, "> expvec.txt" ;
my %G_id ;
my %EXP ;
<I> ;
while(<I>){
  chomp;
  my ($g, $v) = split(/\t/, $_) ;
  my $cnt     = scalar(keys %G_id) ;
  #next if (!(grep {$v eq $_} (@sig_clusts))) ;
  $G_id{$g}   = $cnt ;
  $EXP{$g}    = $v ;
  print O $g, "\t", $v, "\n" ;
}
close I;
print "Done\n";

my @GENES = sort keys %G_id ;

print "Got " . scalar(@GENES) . " genes.\n";

my @MOTIFS = sort keys %M_id ;
my $a_ta_M ;
for (my $i=1 ; $i<=@GENES ; $i++)
{
    $a_ta_M->[$i][0] = $GENES[$i-1] ;
}
for (my $i=1 ; $i<=@MOTIFS ; $i++)
{
    print $MOTIFS[$i-1], "\n" ;
    $a_ta_M->[0][$i] = $MOTIFS[$i-1] ;
}

my $expvec ;
for (my $i=1 ; $i<=@GENES ; $i++)
{
  for (my $j=1 ; $j<=@MOTIFS ; $j++)
    {
      if (defined($status{$MOTIFS[$j-1]}{$GENES[$i-1]}{m}) && $status{$MOTIFS[$j-1]}{$GENES[$i-1]}{m} == 1)
	{
	  $a_ta_M->[$i][$j] = 1 ;
	}
      else
	{
	  $a_ta_M->[$i][$j] = 0 ;
	}
    }
  $expvec->[$i-1] = $EXP{$GENES[$i-1]} ;
}


Subs::saveTable($a_ta_M, "./binarized_table.txt") ;


#
# 
#
my $bayes = Bayes->new ;
$bayes->init($a_ta_M, $expvec, \%status, \%pos_bias, \%ori_bias) ;

my $ta = Table->new ;
$ta->loadFile("./binarized_table.txt") ;
$a_ta_M = $ta->getArray() ;
my $a_ta_H = shift @$a_ta_M ;
shift @$a_ta_H ;

# gene list
my $a_ta_R ;
for (my $i=0 ; $i<@$a_ta_M ; $i++){
    $a_ta_R->[$i] = shift(@{$a_ta_M->[$i]}) ;
}

my $table ;
my $rows ;
for (my $i=0 ; $i<scalar(@$a_ta_R) ; $i++){
  my $g = $a_ta_R->[$i] ;
  $rows->[$i] = $g ;
  for (my $j=0 ; $j<@$a_ta_H ; $j++){
    my $m = $a_ta_H->[$j] ;
    my @v ;
    my $cnt=0 ;
    foreach my $p (@{$status{$m}{$g}{p}}){
      $v[$cnt] = $p*$status{$m}{$g}{o}->[$cnt] ;
    }
    
    #print join("\t", @v) . "\n";
    
    #$table->[$i][$j] = \@v ;
  }
}

#my $result = $bayes->run_bayes_with_biases($rows, $table) ;
my $result = $bayes->run_bayes($rows, $table) ;


open O, "> ./test2.txt" ;
for (my $i=0 ; $i<@$a_ta_R ; $i++){
    print O $a_ta_R->[$i], "\t", $result->[$i] . "\t" . $EXP{$a_ta_R->[$i]} . "\n" ;
}

open O, "> ./test3.txt" ;
my @clust ;
foreach my $g (sort keys %EXP){
  if (! (grep {$_ eq $EXP{$g}} (@clust))){
    push (@clust, $EXP{$g}) ;
  }
}

foreach my $c(@clust){
    my $N1; my $N2; my $N12 ;
    for (my $i=0 ; $i<@$result ; $i++){
	my $g = $a_ta_R->[$i] ;
	if ($EXP{$g} <= $c+2 and $EXP{$g} >= $c-2){
	    $N1++ ;
	}
	if ($result->[$i] <= $c+2 and $result->[$i] >= $c-2){
	    $N2++ ;
	}
	if ($EXP{$g} <= $c+2 and $EXP{$g} >= $c-2 and $result->[$i] <= $c+2 and $result->[$i] >= $c-2){
	    $N12++ ;
	}
    }
    print O $c, "\t", $N1, "\t", $N2, "\t", $N12, "\n" ;
}
