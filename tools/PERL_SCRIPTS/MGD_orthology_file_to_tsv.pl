open IN, $ARGV[0];

my @pos = ();

while (my $l = <IN>) {

  if ($l =~ /^Mouse/) {

    my @s = ("Mouse MGI Acc ID", 
	     "Mouse Chr",
	     "Mouse cM",  
	     "Mouse EntrezGene ID",
	     "Mouse Symbol",
	     "Human Chr", 
	     "Human EntrezGene ID",
	     "Human Symbol",
	     "Mouse Name");

    for (my $i=0; $i<@s; $i++) {
      push @pos, index($l, $s[$i]);
    }
    

  }

  if ($l =~ /^MGI/) {
    
    my $n1 = substr($l, $pos[4], $pos[5]-$pos[4]);
    my $n2 = substr($l, $pos[7], $pos[8]-$pos[7]);
    
    $n1 =~ s/ //g;
    $n2 =~ s/ //g;
    
    print "$n1\t$n2\n";

  }
}
close IN;
