use lib qw(/home/olly/PERL_MODULES);
use Sets;

my $i1 = $ARGV[0];
my $i2 = $ARGV[1];

my %AVG = ();
my %STD = ();
my %POP = ();

for (my $i=$i1; $i<=$i2; $i++) {

    open IN, "7mers_all.txt.$i";

    while (my $l = <IN>) {
	chomp $l;
	my @a = split /\t/, $l, -1;

	push @{ $POP{ $a[0] } }, $a[3]; 
	
	
    }

    close IN;

}


foreach my $k (keys(%POP)) {
    $AVG{$k} = Sets::average($POP{$k});
    $STD{$k} = Sets::stddev($POP{$k});    
}


my %Z = ();
open IN, $ARGV[2];
while (my $l = <IN>) {
    chomp $l;
    my @a = split /\t/, $l, -1;
    #push @{ $POP{ $a[0] } }, $a[3]; 
 
    if ($STD{ $a[0] } > 0.0) {
	my $z = ($a[3] - $AVG{ $a[0] } ) / $STD{ $a[0] };
	print "$a[0]\t$z\n";
    }
}
close IN;

