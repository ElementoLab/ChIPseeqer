use lib qw(/home/olly/PERL_MODULES);
use Sets;

while (my $l = <STDIN>) {
    chomp $l;
    

    
    my @a = split /\t/, $l, -1;

    $a[1] =~ s/_ENST\d+//g;

    my @b = split /\//, $a[1];

    my $a_ref_b = Sets::removeDuplicates(\@b);

    print "$a[0]\t"; print join("/", @$a_ref_b); print "\n";
}
