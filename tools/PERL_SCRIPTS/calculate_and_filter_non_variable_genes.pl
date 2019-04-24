use lib qw(/home/olly/PERL_MODULES);
use Sets;

open IN, $ARGV[0];
my $l = <IN>;



my @TO = ();

my $i = 0;
while (my $l = <IN>) {
    chomp $l;

    my @a = split /\t/, $l;

    my $n = shift @a;

    my $s = Sets::stddev (\@a);
    my $m = Sets::average(\@a);
    
    if (abs($m) < 0.0000001) {
	if ($m < 0.0) {
	    $m = -0.0000001;
	} else {
	    $m =  0.0000001;
	}
    }
    
    my $r = $s; #/$m;
    
    #print "$n\t$r\t$s\t$m\n";
    my @a_tmp = ($i, $r);
    push @TO, \@a_tmp;
    $i++;		 
}
close IN;

@TO = sort { $b->[1] <=> $a->[1] } @TO;

my %H = (); my %H1 = ();
for (my $i=0; $i<$ARGV[1]; $i++) {
    $H { $TO[$i]->[0] } = 1;
    $H1{ $TO[$i]->[0] } = $TO[$i]->[1];
}

open IN, $ARGV[0];
my $l = <IN>;
print $l;
my $i = 0;
while (my $l = <IN>) {
    
    if (defined($H{ $i })) {
	print $l;
    }
    $i++;
}
close IN;
