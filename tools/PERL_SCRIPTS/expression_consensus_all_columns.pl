open IN, $ARGV[0];


my $l = <IN>; chomp $l;
my @a = split /\t/, $l;
my $p = shift @a;
print "$p\t" . join("/", @a) . "\n";

while (my $l = <IN>) {

  
  chomp $l;
  
  my @a = split /\t/, $l, -1;
  my $p = shift @a;
  
  if ($p eq "ID_REF") {
    print "ID_REF\tlogratio\n";
    next; 
  }
  

  my $sum = 0;
  foreach my $r (@a) {
    if (($r ne "nan") && ($r ne "")) {
      $sum += $r;
    }
  }
  
  if ($sum == @a) {
    print "$p\t1\n";    
  } else {
    print "$p\t0\n";
  }
}
