open IN, $ARGV[0];

while (my $l = <IN>) {

        
	chomp $l;
	
	my @a = split /\t/, $l, -1;
	my $p = shift @a;

        if ($p eq "ID_REF") {
           print "ID_REF\tlogratio\n";
           next; 
        }

	my $sum1 = 0;
	for (my $i=0; $i<$ARGV[1]; $i++) {
		$sum1 += $sum1 + $a[$i];
	}
	$sum1 = $sum1 / $ARGV[1];

        my $sum2 = 0;
        for (my $i=$ARGV[1]; $i<$ARGV[1]+$ARGV[2]; $i++) {
                $sum2 += $sum2 + $a[$i];
        }
        $sum2 = $sum2 / $ARGV[2];

	print "$p\t" . sprintf("%3.2f", log($sum1 / $sum2) / log(2.0)) . "\n";
}
