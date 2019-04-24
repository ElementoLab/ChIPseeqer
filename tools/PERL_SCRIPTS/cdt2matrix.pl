my $l = <STDIN>;
my @a = split /\t/, $l, -1;
shift @a;
shift @a;
shift @a;
print "\t" . join("\t", @a);
my $l = <STDIN>;
while (my $l = <STDIN>) {    
    my @a = split /\t/, $l, -1;
    my $n = shift @a;
    shift @a;
    shift @a;
    print $n . "\t" . join("\t", @a);
}
