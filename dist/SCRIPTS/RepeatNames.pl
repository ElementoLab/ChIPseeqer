my %H = ();
open IN, $ARGV[0] or die "Cannot open file";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  
  my @b = split /\-/, $a[0];
  pop @b;
  my $txt =  join("-", @b); 
  if (!defined($H{$txt})) {
    print "$txt\n";
    $H{$txt} = 1;
  }


}
close IN;

