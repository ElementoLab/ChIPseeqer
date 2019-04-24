my $pat = $ARGV[0];
@t = ();
while ($l = <STDIN>) {
    chomp $l;
	
    my @a = split /\t/, $l;
    
    
    foreach my $item (@a) {	
	print "$item\n" if ($item =~ /$pat/);
	
    }
    
    

	
	
	
}




