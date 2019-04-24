use lib qw(/home/olly/PERL_MODULES);

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $h_ref = $ta->getIndexKV(0,1);

my $i = 0;
my @a = ();
while (my $l = <STDIN>) {
    
    next if ($l =~ /^\#/);
    chomp $l;
    
    if ($l =~ /^\>/) {
	($chr, $sta, $end, $fac, $orf) = $l =~ /\>(.+?)\:(\d+)\-(\d+)\|(.+)\-\>(.+)\|PMID/;
    } else {
	my @a_tmp = ($chr, $sta, $end, $fac, $orf, $h_ref->{$orf}, $l);

	push @a, \@a_tmp;
	$i++;
    }
    
    
}

foreach my $r (@a) {
    print join("\t", @$r); print "\n";
}
