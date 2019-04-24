$fin = 0;
$deb = 0;
while (my $l=<STDIN>) {
    
    if ($l =~ /Alignment \(FASTA format\)\:/) {
	$deb = 1; $l=<STDIN>; $l=<STDIN>; $l=<STDIN>;
    }

    if ($l =~ /Sequence tree\:/) {
	$fin = 1; 
    }

    next if ($deb == 0);
    next if ($fin == 1);
    chomp $l;
    $l =~ s/\ //g;
    
    print "$l\n" if ($l ne "");

}
