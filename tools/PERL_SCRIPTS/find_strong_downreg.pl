use lib "$ENV{HOME}/PERL_MODULES";

use Sets;
use strict;

use Getopt::Long;
if (@ARGV == 0) {
  die "Args: --data=FILE --logfold=FLT\n";
}
my $data    = undef;
my $logfold = undef;

GetOptions("data=s"    => \$data,
           "logfold=s" => \$logfold);


open IN, $data;
my $l = <IN>;
#if ($ARGV[1] eq "") {
#print "RefSeq\tCB/NB log(fold-down)\tLY1/NB log(fold-down)\n";
#} else {
  print "RefSeq\tDownCBLY1_NB\n";
#}
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
 #print "$l\n"; 

   if ($a[1] < 1) {
      $a[1] = 1;
    }
    if ($a[2] < 1) {
      $a[2] = 1;
    }
    if ($a[3] < 1) {
      $a[3] = 1;
    }
    my $l1 = Sets::log2($a[2]/$a[1]);
    my $l2 = Sets::log2($a[3]/$a[1]);
 
  my $out = 0;
  if (($logfold <= 0) && (($l1 < $logfold) && ($l2 < $logfold))) {
    $out = 1;
  } elsif (($logfold > 0) && (($l1 > $logfold) && ($l2 > $logfold))) {
    $out = 1;
  } else {
    $out = 0;
  }
  
  print "$a[0]\t$out\n";

}
close IN;

