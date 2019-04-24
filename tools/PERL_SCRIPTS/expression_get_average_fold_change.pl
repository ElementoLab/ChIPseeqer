#!/usr/bin/perl

open IN, $ARGV[0];
my $l = <IN>;

chomp $l;
my @a  = split /\t/, $l, -1;
my $n = shift @a;

my @s1 = ();
for (my $i=0; $i<$ARGV[1]; $i++) {
  push @s1, $a[$i];
}
my $t1 = join("-", @s1);


my @s2 = ();
for (my $i=$ARGV[1]; $i<$ARGV[1]+$ARGV[2]; $i++) {
  push @s2, $a[$i];
}
my $t2 = join("-", @s2);


print "$n\t"; print "$t1/$t2\n";



while (my $l = <IN>) {

        
  chomp $l;
  
  my @a = split /\t/, $l, -1;
  my $p = shift @a;

  if ($p eq "ID_REF") {
    print "ID_REF\tlogratio\n";
    next; 
  }
  
  my $sum1 = 0;
  for (my $i=0; $i<$ARGV[1]; $i++) {
    $sum1 += $sum1 + $a[$i];
  }
  $sum1 = $sum1 / $ARGV[1];
  
  my $sum2 = 0;
  for (my $i=$ARGV[1]; $i<$ARGV[1]+$ARGV[2]; $i++) {
    $sum2 += $sum2 + $a[$i];
  }
  $sum2 = $sum2 / $ARGV[2];
  
  

  #print "$sum1/$sum2\n";

  if ($sum1 == 0.0) {
    $sum1 = 0.00001;
  }  
  if ($sum2 == 0.0) {
    $sum2 = 0.00001;
  }


  my $ratio    = $sum1 / $sum2;

  if ($ARGV[3] ne "") {
    $ratio = 1/$ratio;
  }	

  my $logratio = log ($ratio) / log(2.0);

  

  print "$p\t" . sprintf("%4.3f", $logratio) . "\n";

}
