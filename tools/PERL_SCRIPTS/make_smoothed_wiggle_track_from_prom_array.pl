#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use Sets;
#use strict;

if (@ARGV == 0) {
  die "Args: suffix desc (will add _532.pair)\n";
}

my $suf = $ARGV[0];

my $prob = "/Users/olivier/PEOPLE/WEIMIN/CORREL_RESCUE/ANNOTATION/probes_mapped_to_hg18.txt";
my $desc = $ARGV[1];

my $ta = Table->new;
$ta->loadFile($prob);
my $h_ref_probes = $ta->getIndexShifted(0);


#
# load 532 data
#
my $f532 = "$suf\_532.pair";
if (! -e $f532) {
  $f532 .= ".txt";
}
open IN, $f532 or die "Cannot open $f532\n";
my $l = <IN>;
$l = <IN>;

my %PROMS532 = ();

while (my $l = <IN>) {
  next if ($l =~ /^\#/);
  chomp $l;
  my @a = split /\t/, $l, -1;
  next if ($a[2] =~ /RANDOM/);
  if ($a[9] eq "") {
    die "At $a[3], no intensity.\n";
  }
  my @probes = ($a[3], $a[9]);
  push @{ $PROMS532{$a[2]} }, \@probes;
}
close IN;

#
# load 635 data
#
my $f635 = "$suf\_635.pair";
if (! -e $f635) {
  $f635 .= ".txt";
}
open IN, $f635 or die "Cannot open $f635\n";;
$l = <IN>;
$l = <IN>;

my %PROMS635 = ();

while (my $l = <IN>) {
  next if ($l =~ /^\#/);
  chomp $l;
  my @a = split /\t/, $l, -1;
  next if ($a[2] =~ /RANDOM/);
  my @probes = ($a[3], $a[9]);
  if ($a[9] eq "") {
    die "At $a[3], no intensity.\n";
  }	
  push @{ $PROMS635{$a[2]} }, \@probes;
}
close IN;


#
# calculate ratio profile
#

my @intensities = ();
my %chr_cprobes     = ();  # combine probes
my @nprobes     = ();

foreach my $p (keys(%PROMS532)) {



  my @a_532 = sort sortProbes @{$PROMS532{$p}};
  my @a_635 = sort sortProbes @{$PROMS635{$p}};
  
  my $n = @a_532;

  my @v = ();
  my %CHR = ();

  for (my $i=0; $i<$n; $i++) {
    #print " $a_532[$i]->[0]\n";
    if ($a_532[$i]->[0] ne $a_635[$i]->[0]) {
      die "Problem with data for $p\n";
    }

    my $l = $a_635[$i]->[1] / $a_532[$i]->[1];
    
    my $probe_info = $h_ref_probes->{$a_532[$i]->[0]}; 
    if ($probe_info->[0] ne "") {
      # store chr for this probe
      $CHR{$probe_info->[0]} = 1;
      # make array with pos, ratio
      my @b = ($probe_info->[1], $l, $a_532[$i]->[0]);
    
      # add probe
      push @v, \@b;
    }
  }
  
  # get chromosome
  my $chr = undef;
  my @kk  = keys(%CHR);
  
  if (@kk == 0) {
    #print STDERR "Skipping ($a_532[0]->[0]\n";
  } elsif (@kk > 1) {
    die "multiple chr => " . join("/", @kk) . "\n";
  } else {
    $chr = shift @kk;
  }

  next if ($chr eq "");
  
  #print "Probe set $p, $chr\n";
  
    
  # add probe set
  push @{$chr_cprobes{$chr}}, \@v;
}


#print "track type=wiggle_0 name=\"$desc\" description=\"$desc\"\n";
print "track type=wiggle_0 name=\"$desc\" description=\"$desc\" visibility=\"full\" maxHeightPixels=\"64:64:11\" smoothingWindow=\"10\" viewLimits=\"0:4\" autoScale=\"on\"\n";

foreach my $chr (keys(%chr_cprobes)) {

  print "variableStep chrom=$chr span=10\n";
  #print "1\t0\n";
  #print "11\t4\n";

  # sort probe sets by first probe
  my @a_sorted_ps = sort { $a->[0]->[0] <=> $b->[0]->[0] } @{ $chr_cprobes{$chr} };

  # iterate thru the probes

  my $prev_st = 0;
  my $prev_en = 0;

  foreach my $ps (@a_sorted_ps) {

    # calc st and end of probe set
    my @set = sort { $a->[0] <=> $b->[0] } @$ps;
    my $st  = $set[0]->[0];
    my $en  = $set[$#set]->[0]+50;

    if (!Sets::sequencesOverlap($st, $en, $prev_st, $prev_en)) {
      &smoothe_ps($ps);
      $prev_st = $st;  
      $prev_en = $en;
    }
    


  }

}


sub smoothe_ps {
  my ($r) = @_;
  
  my @set = sort { $a->[0] <=> $b->[0] } @$r;

  
  my $st  = $set[0]->[0];
  my $en  = $set[$#set]->[0]+50;
  my $off = $st;
  
  print STDERR "$st\t$en\t";
  print STDERR scalar(@set) . "\t" . $set[0]->[2] . "\n";
  
  if ( ($en - $st) > 10000) {
    return;
  }


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
    print "$ii\t" . sprintf("%4.3f", $a_newint[$i]) . "\n";
  

  }
}



sub sortProbes {
  my ($p1) = $a->[0] =~ /^CHR[XY\d]+PR{0,1}0*(\d+?)$/;
  my ($p2) = $b->[0] =~ /^CHR[XY\d]+PR{0,1}0*(\d+?)$/;
  #print "$p1 <=> $p2\n";
  return $p1 <=> $p2;
}
