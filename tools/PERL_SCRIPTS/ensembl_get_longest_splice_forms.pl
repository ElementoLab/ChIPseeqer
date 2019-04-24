open IN, $ARGV[0];

my %MRNA_LEN = ();

while (my $l = <IN>) {
    chomp $l;
    
    my @a = split /\t/, $l, -1;

    $NAME   { $a[7] } = 1;
    $CHR    { $a[7] } = $a[1];
    $STRAND { $a[7] } = $a[4]; 
    
    # MRNA
    if (defined($M_STA{ $a[7] })) {

      # old mRNA length
      my $old_l = $M_END{ $a[7] } - $M_STA{ $a[7] };
      
      # new mRNA length
      my $new_l = $a[6]           - $a[5];
      
      if ($new_l > $old_l) {
	# new mRNA larger, update mRNA coordinates
	$M_STA{ $a[7] } = $a[5];
	$M_END{ $a[7] } = $a[6];

	# also update coding sequence coordinates	
	$E_STA{ $a[7] } = $a[2];
	$E_END{ $a[7] } = $a[3];
      }
            
    } else {
      
      # create mRNA coords      
      $M_STA{ $a[7] } = $a[5];
      $M_END{ $a[7] } = $a[6];
      
      # create corresp coding coords
      $E_STA{ $a[7] } = $a[2];
      $E_END{ $a[7] } = $a[3];
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
