while (my $l = <STDIN>) {
    chomp $l;
    my @a = split /\t/, $l, -1;
    
    $a[0] =  lc($a[0]);
    $a[0] =~ s/_//g;
    
    print join("\t", @a); print "\n";
}
