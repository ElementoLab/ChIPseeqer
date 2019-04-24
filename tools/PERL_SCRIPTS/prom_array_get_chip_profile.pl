#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use Table;
use Sets;

use strict;

my $hsap = $ARGV[0];
my $suff = $ARGV[1];
my $prob = "/Users/olivier/PEOPLE/WEIMIN/CORREL_RESCUE/ANNOTATION/probes_mapped_to_hg18.txt";
my $desc = $ARGV[2];

my $ta = Table->new;
$ta->loadFile($prob);
my $h_ref_probes = $ta->getIndexShifted(0);


my $a_ref1 = &get_probes_info($hsap, "$suff\_532.pair");
my $a_ref2 = &get_probes_info($hsap, "$suff\_635.pair");

my %H = ();
foreach my $r (@$a_ref2) {
  $H{ $r->[0] } = $r->[1];
  #print join("\t", @$r) . "\n";
}
#print "\n";

my @set = ();
my %CHR = ();
foreach my $r (@$a_ref1) {

  # 635/532
  # my $l = log($H{$r->[0]}/$r->[1]);
  my $l = $H{$r->[0]}/$r->[1];

  my $probe_co = $h_ref_probes->{$r->[0]}; 
 # print  $probe_co->[1] . "\t$l\n";  

  $CHR{$probe_co->[0]} = 1;
  my @b = ($probe_co->[1], $l);
  push @set, \@b;
}

@set = sort { $a->[0] <=> $b->[0] } @set;


my $st  = $set[0]->[0];
my $en  = $set[$#set]->[0]+50;
my $off = $st;

#print STDERR "$st\t$en\t";
#print STDERR scalar(@set) . "\n";

# fill up array
my @a_int = ();
my @a_num = ();
foreach my $s (@set) {
  for (my $i=$s->[0]; $i<$s->[0]+50; $i++) {
    $a_int[$i-$off] += $s->[1];
    $a_num[$i-$off] ++;
  }
}

# average
for (my $i=0; $i<$en-$off; $i++) {
  if ($a_num[$i] > 1) {
    $a_int[$i] /= $a_num[$i];
  }
}

my $chr = undef;
my @kk  = keys(%CHR);

if (@kk > 1) {
  die "multiple chr.\n";
} else {
  $chr = shift @kk;
}

print "track type=wiggle_0 name=\"$desc\" description=\"$desc\"\n";
print "variableStep chrom=$chr span=10\n";
print "1\t0\n";
print "11\t4\n";
my @postval = ();
my $prevpostval = undef;
my $prevpostpos = undef;
for (my $i=$en-$off-1; $i>=0; $i--) {
  if ($a_int[$i] ne "") {
    $prevpostval = $a_int[$i];
    $prevpostpos = $i;
  }
  $postval[$i] = [ $prevpostval, $prevpostpos ];
}

my @a_newint = ();
my $preval = undef;
my $prepos = undef;
for (my $i=0; $i<$en-$off; $i+=10) {

  if ($a_int[$i] ne "") {
    $a_newint[$i] = $a_int[$i];
    $preval = $a_int[$i];
    $prepos = $i;
  } else {
    
    my $a = ($postval[$i]->[0] - $preval) / ($postval[$i]->[1] - $prepos);
    my $b = $preval - $a * $prepos;

    $a_newint[$i] = $a * $i + $b;

  }

  my $ii = $i + $off;
  print "$ii\t$a_newint[$i]\n";
  

}



sub get_probes_info {
  my ($id, $suff) = @_;

  open IN, "$suff" or die "Cannot open $suff\n";
  
  my @c = ();
  while (my $l = <IN>) {
    chomp $l;
    my @a = split /\t/, $l, -1;
    
    if ($a[2] eq $hsap) {
      
      my @b = ($a[3],$a[9]);
      
      push @c, \@b;
      #print join("\t", @b) . "\n";

    }

  }
  
  close IN;

  return \@c;
}
