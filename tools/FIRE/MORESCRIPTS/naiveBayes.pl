use lib "$ENV{FIREDIR}/SCRIPTS";
use Getopt::Long;
use Sets ;
use Data::Dumper ;

use Table ;
use strict;
use Bayes ;

my $expfile  = undef ;
my $profiles = undef ;
my $testfile = undef;

GetOptions ('expfile=s'     => \$expfile,
	    'testfile=s'    => \$testfile);


my $expfile_f   = Sets::filename($expfile);
my $namesfile   = "$expfile\_FIRE/DNA/$expfile_f.motifnames";
my $summaryfile = "$expfile\_FIRE/DNA/$expfile_f.summary";
my $distfile    = "$expfile\_FIRE/DNA/$expfile_f.dist.rep";
my $profiles    = "$expfile\_FIRE/DNA/$expfile_f.profiles";


#
# reading profiles
#
my %status ;
open I, "< $profiles" or die "Couldn't open $profiles";
while(<I>){
    chomp;
    my ($m, $g, $p, $o, @a) = split(/\t/, $_);
    $status{$m}{$g}{m} = 1 ;
    push(@{$status{$m}{$g}{p}}, $p) ;
    push(@{$status{$m}{$g}{o}}, $o) ;
}
close I;


open I, "< $summaryfile" or die "Couldn't: $summaryfile";
my @sig_clusts ;
my %pos_bias ;
my %ori_bias ;
my %M_id;
my $cnt = 0;
while(<I>){
  chomp;
  my @a = split(/\t/, $_) ;
  $M_id    {$a[0]} = $cnt ;
  $pos_bias{$a[0]} = $a[9] ;
  $ori_bias{$a[0]} = $a[10] ;
  for (my $i=12 ; $i<@a ; $i++){
    if (!(grep {$_ eq $a[$i]} (@sig_clusts))){
      push(@sig_clusts, $a[$i]) ;
    }
  }
  $cnt++;
}
close I;

open I, "< $distfile" or die "Couldn't: $distfile";
while(<I>){
    chomp;
    my ($m, $l, $e, @a) = split(/\t/, $_) ;
    if ($l eq "d_avg" && $e eq "nan"){
	pop @a ;
	$status{$m}{bins} = \@a ;
    }
}
close I;


my @MOTIFS    = sort keys %M_id ;
my $n_motifs  = @MOTIFS;

my $a_ref_res = makeFeatureMatrix($expfile, \%status);

my $a_ref_tes = makeFeatureMatrix($testfile, \%status);

#
# Bayes 
#
my $bayes = Bayes->new ;

$bayes->setFeatureNames (\@MOTIFS);
$bayes->setFeatureMatrix($a_ref_res->[0]);
$bayes->setClassVector  ($a_ref_res->[1]);
$bayes->setObjectNames  ($a_ref_res->[2]);

$bayes->train();  # learn probabilities

my $classes = $bayes->classify($a_ref_tes->[0]);  # classify (using the same features as the one used to learn)

my $n_genes = scalar( @{$a_ref_tes->[0]} );
print "Gene\tREAL\tPRED\n";
for (my $i=0; $i<$n_genes; $i++) {
  print "$a_ref_tes->[2]->[$i]\t$a_ref_tes->[1]->[$i]\t$classes->[$i]\n";
}




sub makeFeatureMatrix {

  my ($expfile, $prof) = @_;

  open I, "< $expfile" or die "Couldn't: $expfile" ;
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
  }
  close I;

  my @GENES  = sort keys %G_id ;
  my $n_genes  = @GENES;

  #print "Got $n_genes genes and $n_motifs motifs.\n";
  
  my $a_ref_C = [];
  for (my $i=0 ; $i<@GENES ; $i++) {
    $a_ref_C->[$i] = $EXP{ $GENES[$i] };
  }
  
  #
  # make a GENES x MOTIFS matrix
  #
  
  my $a_ref_F = [];
  
  for (my $i=0 ; $i<@GENES ; $i++) {
    for (my $j=0 ; $j<@MOTIFS ; $j++) {
      if (defined($prof->{$MOTIFS[$j]}{$GENES[$i]}{m}) && ($prof->{$MOTIFS[$j]}{$GENES[$i]}{m} == 1)) {
	$a_ref_F->[$i][$j] = 1 ;
      } else {
	$a_ref_F->[$i][$j] = 0 ;
      }
    }
  }

  return [ $a_ref_F, $a_ref_C, \@GENES ];

}
