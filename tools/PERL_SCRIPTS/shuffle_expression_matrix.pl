use lib qw(/home/olly/PERL_MODULES);
use Sets;

srand;


open IN, $ARGV[0];
my $l = <IN>;
print $l;
while (my $l = <IN>) {
    chomp $l;
    my @a = split /\t/, $l, -1;
    my $n = shift @a;
    my $a_ref_shu = Sets::shuffle_array(\@a);
    print "$n\t" . join("\t", @$a_ref_shu) . "\n";
}
close IN;



