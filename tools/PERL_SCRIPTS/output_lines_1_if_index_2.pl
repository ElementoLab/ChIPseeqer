BEGIN{ $home = `echo \$HOME`; chomp $home}

use lib "$home/PERL_MODULES";

use Sets;

my $h_ref = Sets::getIndex($ARGV[1]);


open IN, $ARGV[0];
my $l = <IN>;
print "$l";
while (my $l = <IN>) {
    chomp $l;

    my @a = split /\t/, $l, -1;
    
    if (defined($h_ref->{ $a[ $ARGV[2] ] })) {
	#$a[0] = $h_ref->{ $a[0] } . "_" . $a[0];
	print join("\t", @a) . "\n";
    }

}
close IN;


