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

    # test what is 
    my $a1 = $a->[$cc];
    my $b1 = $b->[$cc];


    #print "Compare $a1 and $b1\n";
    
    # test whether $a1 and $b1 are integer 
    my $a1i = 0;
    my $b1i = 0;

    
    $a1i = 1 if ($a1 !~ /[A-Za-z\_]/);
    $b1i = 1 if ($b1 !~ /[A-Za-z\_]/);
    #print sprintf("%d", $a1); print "\t";
    
    #print "$a1i .. $b1i ..\n";

    if (($a1i == 1) && ($b1i == 1)) { 	
	return $a1 <=> $b1;
    }
    
    elsif ( ($a1i != 1) && ($b1i != 1) ) { 	
	return $a1 cmp $b1;
    }
    
    elsif ( ($a1i == 1) && ($b1i == 0) ) {
	return -1;
    }

    else {
	return 1;
	
    }
}



