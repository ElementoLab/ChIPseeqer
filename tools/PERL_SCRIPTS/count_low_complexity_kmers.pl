use lib qw(/home/olly/PERL_MODULES);
use Sets;
use Table;

#Sets::printSet(Sets::allkmers($ARGV[0]));

#Sets::printSet(

my $ta = Table->new;
$ta->loadFile($ARGV[0]);


my $a_ref = $ta->getArray();

my %h = ();
my $C = \%h;
my $cnt = undef;
my $l = 0;
foreach my $r (@$a_ref) {

    my $H = Sets::countNucleotides($r->[0]);

    $l = length($r->[0]);

    if ($H->{"A"} >= length($r->[0])-2) { 
	$C->{"A"}->[$H->{"A"}] += $r->[1];
	$cnt->{"A"}->[$H->{"A"}]++;
    }

    if ($H->{"T"} >= length($r->[0])-2) { 
	$C->{"A"}->[$H->{"T"}] += $r->[1];
	$cnt->{"A"}->[$H->{"T"}]++;

    }
    
    if ($H->{"C"} >= length($r->[0])-2) { 
	$C->{"C"}->[$H->{"C"}] += $r->[1];	
	$cnt->{"C"}->[$H->{"C"}]++;
    }

    if ($H->{"G"} >= length($r->[0])-2) { 
	$C->{"C"}->[$H->{"G"}] += $r->[1];
	$cnt->{"C"}->[$H->{"G"}]++;	
    }
    
}

foreach my $k (keys(%$C)) {
    my $i = -1;
    foreach my $r (@{ $C->{$k} }) {
	$i++;
	next if (!defined($r));
	my $o = int(0.5+$r/$cnt->{$k}->[$i]);
	
	print "$k $i/$l\t$o\n";
	
    }
}


