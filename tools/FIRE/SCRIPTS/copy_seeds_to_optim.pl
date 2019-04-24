use lib "$ENV{FIREDIR}/SCRIPTS";

&copySeedsToOptim($ARGV[0], $ARGV[1]);


#
# copy seed format to optim format
#
sub copySeedsToOptim {
  my ($seeds7, $optim7) = @_;
  
  open OOU, ">$optim7" or die "copySeedsToOptim: Cannot open $optim7 for writing.\n";
  open OOO, $seeds7 or die "copySeedsToOptim: Cannot open $seeds7.\n";
  while (my $l = <OOO>) {
    chomp $l;
    my @a = split /\t/, $l, -1;
    print OOU "$a[0]\t$a[1]\t0\t0.0\t$a[4]\t$a[5]\n";
  }
  close OOO;
  close OOU;
  
}
