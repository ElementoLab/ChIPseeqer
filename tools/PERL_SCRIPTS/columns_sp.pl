
@t = ();
while ($l = <STDIN>) {
    chomp $l;
    my @a = split /[\t\ ]/, $l;
	my @l = ();
	foreach my $c (@ARGV) {
		push @l, $a[$c-1];
		
	}
	print join ("\t", @l) . "\n";
}




