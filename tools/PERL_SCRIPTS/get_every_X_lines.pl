my $i = 0;
while (my $l = <STDIN>) {
    
    if ($i % 10 == 0) {
	print $l;
    }

    $i ++;
}
