open OUT, "t";

while ($l = <STDIN>) {
    
    chomp $l;

    if ($l =~ /Cluster/) {
	
	close OUT;
	
	$l =~ s/\ /_/g;

	open OUT, ">$l.txt";
    } else {
	
	print OUT "$l\n";
    }
}


close OUT;
