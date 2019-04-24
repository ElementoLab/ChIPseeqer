srand;
    @new = ();
    @old = 1 .. 10;  # just a demo
    for( @old ){
        my $r = rand @new+1; print "r=$r, new[r]=" . $new[$r] . "\n"; 
        push(@new,$new[$r]);
        $new[$r] = $_; 

        print join("\t", @new); print "\n";
    }
