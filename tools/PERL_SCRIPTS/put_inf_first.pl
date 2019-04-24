my $s = "";
my $i = "";
while (my $l = <STDIN>) {
   
    if ($l =~ /inf/) {
	$i .= $l;
    } else {
	$s .= $l;
    }
    
    
}

print "$i$s";
