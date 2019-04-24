my %H = ();
while (my $l = <STDIN>) {
    
    chomp $l;
    
    my @a = split /\t/, $l, -1;
    
    $a[0] =~ s/\-P.+$//g;
    
    if (defined($H{ $a[0] })) {
	my $s1 = abs($a[3] - $a[2]);

	my $s2 = abs($H{$a[0]}->[3] - $H{$a[0]}->[2]);
	
	if ($s1 > $s2) {
	    $H{ $a[0] } = \@a;
	}
    } else {
	$H{ $a[0] } = \@a;
    }
    
    
    
    
}


reset(%H);
while (my ($k, $v) = each(%H)) {
    print join("\t", @$v) . "\n"; 
}
