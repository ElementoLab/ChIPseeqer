my @a = <STDIN>;

$yes = 1;
foreach $aa (@a) {
	
	if ($yes == 1) {
		print $aa;
		$yes = 0;
	} else {
		$yes = 1;
	}

}

