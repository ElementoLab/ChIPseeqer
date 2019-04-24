#!/usr/bin/perl
$cc = $ARGV[0];

@t = ();
while ($l = <STDIN>) {
    chomp $l;
    my @a = split /\t/, $l, -1;
    push @t, \@a;
}



foreach $r (sort sortLines @t) {
    print join("\t", @$r) . "\n";
}


sub sortLines {

    return $a->[$cc] <=> $b->[$cc];

}



