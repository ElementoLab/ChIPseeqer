#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;
use Getopt::Long;
use strict;

my $expt = undef;
my $cont = undef;
my $doavg = 1;
my $dolog = 0;

GetOptions("expt=s" => \$expt,
           "cont=s" => \$cont,
	   "doavg=s"  => \$doavg,
	   "dolog=s"  => \$dolog);


my $nmcol = 7;
my $sicol = 8;
my $ncski = 9;

my @a_all_meds = ();
my @a_cont     = ();    

if (defined($cont)) {
  # control
  
  # open file
  open IN, $cont or die "Cannot open $cont\n";
  for (my $i=0; $i<$ncski; $i++) {
    my $l = <IN>;
  }

  my $n = 0;
  my @a_v = ();
  while (my $l = <IN>) {
    chomp $l;
    my @a = split /\t/, $l;
    
    my $NM = $a[$nmcol];
    next if ($NM !~ /^NM/);
    my $VA = $a[$sicol];
    
    #print "$NM\t$VA\n";
    my @b = ($NM, $VA);
    push @a_cont, \@b;
    push @a_v, $VA;
    $n++;
    
  }
  close IN;
  
  my $medcont = Sets::median(\@a_v);
  #print "# CON = $medcont\n";
  
  push @a_all_meds, $medcont;
}



my @files = split /\,/, $expt;

my @a_expts = ();


foreach my $f (@files) {
  
  # open file
  open IN, $f or die "Cannot open $f\n";
  my @a_v = ();

  for (my $i=0; $i<$ncski; $i++) {
    my $l = <IN>;
  }
  my @a_e;
  my $n = 0;
  while (my $l = <IN>) {
    chomp $l;
    my @a = split /\t/, $l;
    
    my $NM = $a[$nmcol];
    next if ($NM !~ /^NM/);
    my $VA = $a[$sicol];

    #print "$NM\t$VA\n";
    my @b = ($NM, $VA);
    push @a_e, \@b;
    push @a_v, $VA;
    $n++;
  }
  close IN;

  my $med = Sets::median(\@a_v);
  push @a_all_meds, $med;
  print "# MED = $med\n";

  push @a_expts, \@a_e;
  
}

my $avg_med = Sets::average(\@a_all_meds);

my @a_scale = @a_all_meds;
foreach my $r (@a_scale) {
  $r = $avg_med / $r;
}


my $nf = @files;
my $ni = @{ $a_expts[0] }; 

print "GENE\tFOLD\n";

my @data = ();


my $sca = undef;
if (defined($cont)) {
  $sca = shift @a_scale;
}

for (my $i=0; $i<$ni; $i++) {
  

  if (defined($cont)) {
    $a_cont[$i]->[1] *= $sca;
    print "$a_cont[$i]->[0]";
    push @{ $data[0] }, $a_cont[$i]->[1];
  } else {
    print  $a_expts[0]->[$i]->[0];  
  }


  my $avgl = 0;
  my $txt  = ""; 
  for (my $j=0; $j<$nf; $j++) {
    #die "Problem with NMs.\n" if ($a_cont[$i]->[0] ne $a_expts[$j]->[$i]->[0]);

    $a_expts[$j]->[$i]->[1] *= $a_scale[$j];

    my $r = $a_expts[$j]->[$i]->[1];

    if (defined($cont)) {
      $r /= $a_cont[$i]->[1];
    }

    $avgl += $r;

    if (defined($cont)) {
      push @{ $data[$j+1] }, $a_expts[$j]->[$i]->[1];      
    } else {
      push @{ $data[$j] }, $a_expts[$j]->[$i]->[1];      
    }
    
    
    if ($dolog == 1) {
      $a_expts[$j]->[$i]->[1] = Sets::log2($a_expts[$j]->[$i]->[1]);
    }

    $txt .= "\t" . $a_expts[$j]->[$i]->[1];

  }
  
  #for (my $j=0; $j<$nf; $j++) {
  
  #}
  $avgl /= $nf;
  
  if ($doavg == 1) {
    print "\t$avgl\n";
  } else {
    print "$txt\n";
  }
  

  
}
  
foreach my $r (@data) {
  print "# " . Sets::median($r) . "\n";
}    
