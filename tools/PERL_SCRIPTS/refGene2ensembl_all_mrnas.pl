open IN, $ARGV[0];

my %MRNA_LEN = ();

while (my $l = <IN>) {

  chomp $l;
    
  my @a = split /\t/, $l, -1;
  
  #$NAME   { $a[1] } = 1;
  #$CHR    { $a[1] } = $a[2];
  my $st = ($a[2] eq "+"?1:-1);
  
  #$M_STA  { $a[1] } = $a[4];
  #$M_END  { $a[1] } = $a[5];
  
  #$E_STA  { $a[1] } = $a[6];
  #$E_END  { $a[1] } = $a[7];
 
  print "$a[0]\t$a[1]\t$a[5]\t$a[6]\t$st\t$a[3]\t$a[4]\n";
 
}

close IN;

