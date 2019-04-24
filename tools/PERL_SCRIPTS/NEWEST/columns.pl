#!/usr/bin/perl
@t = ();
while ($l = <STDIN>) {
    chomp $l;
    my @a = split /\t/, $l, -1;
	my @l = ();
	foreach my $c (@ARGV) {
		push @l, $a[$c];
		
	}
	print join ("\t", @l) . "\n";
}




