use lib qw(/home/olly/PERL_MODULES);

open IN, $ARGV[0];
my $l = <IN>;
my @a = split /\t/, $l, -1;
shift @a; shift @a; shift @a; 
print "\t"; print join("\t", @a); 
$l = <IN>;
while (my $l = <IN>) {
    chomp $l;
    my @a = split /\t/, $l, -1;
    shift @a;
    
    my $g = shift @a; shift @a;
    my @b = split /\ \|\| /, $g, -1;

    if ($b[1]) {
	print "$b[1]\t";
	print join("\t", @a); print "\n";
    }
}
close IN;
