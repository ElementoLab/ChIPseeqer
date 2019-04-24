open IN, $ARGV[0];

my %MRNA_LEN = ();

while (my $l = <IN>) {
    chomp $l;
    
    my @a = split /\t/, $l, -1;

    $NAME   { $a[1] } = 1;
    $CHR    { $a[1] } = $a[2];
    $STRAND { $a[1] } = ($a[3] eq "+"?1:-1);
    
    # MRNA
    if (defined($M_STA{ $a[1] })) {

      # old mRNA length
      my $old_l = $M_END{ $a[1] } - $M_STA{ $a[1] };
      
      # new length
      my $new_l = $a[5]           - $a[4];
      
      if ($new_l > $old_l) {
	# new mRNA larger, update mRNA coordinates
	$M_STA{ $a[1] } = $a[4];
	$M_END{ $a[1] } = $a[5];

	# also update coding sequence coordinates	
	$E_STA{ $a[1] } = $a[6];
	$E_END{ $a[1] } = $a[7];
      }
            
    } else {
      
      # create mRNA coords      
      $M_STA{ $a[1] } = $a[4];
      $M_END{ $a[1] } = $a[5];
      
      # create corresp coding coords
      $E_STA{ $a[1] } = $a[6];
      $E_END{ $a[1] } = $a[7];
    }
    
}

close IN;

foreach my $g (keys(%NAME)) {
  print "$g\t";
  print "$CHR{$g}\t";
  print "$E_STA{$g}\t";
  print "$E_END{$g}\t";
  print "$STRAND{$g}\t";
  print "$M_STA{$g}\t";
  print "$M_END{$g}\n";
}
