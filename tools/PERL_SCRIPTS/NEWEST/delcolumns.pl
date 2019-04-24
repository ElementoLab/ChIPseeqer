#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
#use lib qw(/home/elemento/PERL_MODULES); 
use lib "$home/PERL_MODULES";

use Sets;

@t = ();
while ($l = <STDIN>) {
    chomp $l;
    my @a = split /\t/, $l, -1;
	my @l = ();
	for (my $i=0; $i<scalar(@a); $i++) {
                $c = $a[$i];
		push @l, $c if (!Sets::in_array($i, @ARGV));
		
	}
	print join ("\t", @l) . "\n";
}




