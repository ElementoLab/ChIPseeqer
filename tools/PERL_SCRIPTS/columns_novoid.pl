#!/usr/bin/perl
@t = ();
while ($l = <STDIN>) {
    chomp $l;
	
    my @a = split /\t/, $l;
	my @l = ();
	my $empty = 0;
	foreach my $c (@ARGV) {

		
		push @l, $a[$c-1];

		my $t = $a[$c-1]; $t =~ s/\s//g;
	
		$empty ++ if ($t eq "");
	}
	

	next if (scalar(@l) == $empty);
	
	print join ("\t", @l) . "\n";
}




