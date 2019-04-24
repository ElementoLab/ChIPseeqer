#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";
use Sets;
use strict;

# read urvival
open IN, $ARGV[2] or die "Cannot open $ARGV[2]\n";
my %SU = ();
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  $SU{ $a[4] } = [ $a[13], $a[16] ];
}
close IN;



my $num = 50;

# read set of probes 

my %PROBES1 = ();
open IN, $ARGV[0] or die "Cannot open $ARGV[0]\n";
my $l = <IN>;
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  $PROBES1{$a[0]} = 1;
}
close IN;


my %PROBES2 = ();
open IN, $ARGV[1] or die "Cannot open $ARGV[1]\n";
my $l = <IN>;
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  $PROBES2{$a[0]} = 1;
}
close IN;

my @aset1 = keys(%PROBES1);
my @aset2 = keys(%PROBES2);
my $a_ref_int = Sets::getOverlapSet(\@aset1, \@aset2);
my %IDXP = ();
print scalar(@$a_ref_int) . " probes\n";
foreach my $i (@$a_ref_int) {
  $IDXP{$i} = 1;
}



# 1. load exp1 and exp2, then get ranks
#get_rank_from_array

my @a_probes = ();
my @a_exp1 = ();
my @a_exp2 = ();

open IN, $ARGV[0] or die "Cannot open $ARGV[0]\n";
my $l = <IN>;
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  
  if (defined($IDXP{$a[0]})) {
    push @a_probes, $a[0];
    push @a_exp1,   $a[1];
    push @a_exp2,   $a[2];
  }

}
close IN;

my $a_exp1_ranks = Sets::get_rank_from_array(\@a_exp1);
my $a_exp2_ranks = Sets::get_rank_from_array(\@a_exp2);

# read table 2
open IN, $ARGV[1] or die "Cannot open $ARGV[1]\n";
my $l = <IN>; chomp $l;
my @SAMPLES = split /\t/, $l;
shift @SAMPLES;

my @COLS = ();
my @GENES = ();
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  my $g = shift @a;
  if (defined($IDXP{$g})) {
    for (my $i=0; $i<@a; $i++) {
      push @{$COLS[$i]}, $a[$i];      
    }
    push @GENES, $g;
  }
}
close IN;

my $ncols = @COLS;
for (my $i=0; $i<$ncols; $i++) {
  next if ($SAMPLES[$i] eq "");
  my $a_exp3_ranks = Sets::get_rank_from_array($COLS[$i]);
  # reorder
  my %h_ranks = ();
  my %h_vals  = ();
  for (my $j=0; $j<@$a_exp3_ranks; $j++) {
    $h_ranks{$GENES[$j]} = $a_exp3_ranks->[$j];
    $h_vals {$GENES[$j]} = $COLS[$i][$j];
  }
  
  my @a_q = ();
  my @a_c = ();
  my @a_k = ();

  for (my $j=0; $j<$num; $j++) {
    #print "$a_probes[$j]\t$a_exp1_ranks->[$j]($a_exp1[$j])\t$a_exp2_ranks->[$j]($a_exp2[$j])\t" . $h_ranks{$a_probes[$j]} . "(" . $h_vals{$a_probes[$j]} . ")\n";    
    push @a_q, $a_exp1_ranks->[$j];
    push @a_c, $a_exp2_ranks->[$j];
    push @a_k, $h_ranks{$a_probes[$j]};

  }
  for (my $j=@a_probes-$num; $j<@a_probes; $j++) {
    #print "$a_probes[$j]\t$a_exp1_ranks->[$j]($a_exp1[$j])\t$a_exp2_ranks->[$j]($a_exp2[$j])\t" . $h_ranks{$a_probes[$j]} . "(" . $h_vals{$a_probes[$j]} . ")\n";    

    push @a_q, $a_exp1_ranks->[$j];
    push @a_c, $a_exp2_ranks->[$j];
    push @a_k, $h_ranks{$a_probes[$j]};
  }

  #my $t = Sets::shuffle_array(\@a_k);
  #@a_k = @$t;

  my $r_qk = Sets::pearson(\@a_q, \@a_k);
  my $r_kc = Sets::pearson(\@a_k, \@a_c);
  my $r_qc = Sets::pearson(\@a_q, \@a_c);

  my $pe   = ( $r_qk - $r_kc * $r_qc) / ( sqrt(1 - $r_kc * $r_kc) *  sqrt(1 - $r_qc * $r_qc) );
  my $sa = $SAMPLES[$i];
  print "$sa\t$pe\t( r_qk=$r_qk, r_kc=$r_kc, r_qc=$r_qc)\t$SU{$sa}->[0]\t$SU{$sa}->[1]\n";
  #last;
}








