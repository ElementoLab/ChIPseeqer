my $l = <STDIN>;
my @a = split /\t/, $l, -1;
shift @a;
$l = "ORF\tDescription\tGWEIGHT\t" . join("\t", @a);
map { $_ = 1 } @a;
print $l;
print "EWEIGHT\t1\t1\t" . join("\t", @a) . "\n";
while (my $l = <STDIN>) {
    chomp $l;
    my @a = split /\t/, $l, -1;
    my $o = shift @a;
    $l  = "$o\tGene\t1\t"; 
    $l .= join("\t", @a);
    print "$l\n";
}
