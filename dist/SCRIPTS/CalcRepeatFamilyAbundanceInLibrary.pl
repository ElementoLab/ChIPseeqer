use strict;
open IN, $ARGV[0] or die "Cannot open file";

my @fam = ();
if ($ARGV[1] ne "") {
  open IN1, "$ENV{CHIPSEEQERDIR}/DATA/RepMask3.2.7_annotation_families.txt" or die "Cannot open file";
  while (my $l = <IN1>) {
    chomp $l;
    my @a = split /\t/, $l, -1;
    next if ($a[0] =~ /\?$/);
    next if ($a[0] =~ /Simple_repeat/);
    next if ($a[0] =~ /Low_complexity/);

    push @fam, $a[0];
  }
  close IN1;
}
#print @fam;
my %C = ();

while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  
  if ($ARGV[1] eq "") {
     #print "shu]n";  
    if ($a[0] =~ /^Alu/) {
      $C{"Alu"} += $a[3];
    }
    
    if ($a[0] =~ /SINE/) {
      $C{"SINE"} += $a[3];
    }
    
    if ($a[0] =~ /LINE/) {
      $C{"LINE"} += $a[3];
    }
    
    if ($a[0] =~ /L1/) {
      $C{"L1"} += $a[3];
    }
    
    if ($a[0] =~ /L2/) {
      $C{"L2"} += $a[3];
    }
    
    
    if ($a[0] =~ /LTR/) {
      $C{"LTR"} += $a[3];
    }
    
    # indic
    $C{$a[0]} = $a[3];
  
  } else {
  #print "wtf?\n";   
    foreach my $t (@fam) {
      # print "does $a[0] match /\-$t$/\n";
      if ($a[0] =~ /\-$t$/) {
	$C{$t} += $a[3];
      }
    }
    
  }
}
close IN;

my @rep = ();
if ($ARGV[1] eq "") {
  @rep = keys(%C);
} else {
  @rep = @fam;
}


my $ff = $ARGV[0];
$ff =~ s/\_combined\_mastertable.+$//;
$ff =~ s/\_master\_table.+$//;

print "REP\t$ff\n";
foreach my $k (@rep) {
  if (!defined($C{$k})) {
    $C{$k} = 0.0;
  }
  $C{$k} = sprintf("%3.2f", $C{$k});
  print "$k\t$C{$k}\n";
}
