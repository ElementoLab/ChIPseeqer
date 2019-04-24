use lib qw(/home/olly/PERL_MODULES);



$cc = $ARGV[0];
$li = $ARGV[1];

@t = ();
while ($l = <STDIN>) {
    chomp $l;
    my @a = split /\t/, $l;

    #$a[3] = $a[1] / (10*Sets::log2($a[2]));

    push @t, \@a if ($a[1] < $li);
}



foreach $r (sort sortLines @t) {
    print join("\t", @$r) . "\n";
}


sub sortLines {

    return $b->[$cc] <=> $a->[$cc];

}



