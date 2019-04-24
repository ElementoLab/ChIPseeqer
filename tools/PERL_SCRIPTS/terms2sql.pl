my @a = ();
while (my $l = <STDIN>) {

    chomp $l;


    next if ($l =~ /\!/);


    my @a = split /\t/, $l, -1;

    $a[0] =~ s/GO\://g;
    $a[0] = int($a[0]);
    
    print "$a[0]\t$a[1]\t\n";

}

    
