#!/usr/bin/perl
use lib qw(/home/elemento/PERL_MODULES);
use Sets;

my $l = <STDIN>;
chomp $l;
my @a = split /\t/, $l;
my @todel = ();
my @toshow = ();
my $i = 0;
foreach my $r (@a) {
    if (Sets::in_array($r, @ARGV)) {
	push @todel, $i; 
    } else {
	push @toshow, $r;
    }
    $i++;
}
print join("\t", @toshow); print "\n";


@t = ();
while ($l = <STDIN>) {
    chomp $l;
    my @a = split /\t/, $l, -1;
    my @l = ();
    for (my $i=0; $i<scalar(@a); $i++) {
	$c = $a[$i];
	push @l, $c if (!Sets::in_array($i, @todel));
	
    }
    print join ("\t", @l) . "\n";
}




