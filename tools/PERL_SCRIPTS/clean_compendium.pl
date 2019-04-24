$l = <STDIN>;
print uc($l);
while ($l = <STDIN>) {
    
    chomp $l;
    
    my @a =  split /\t/, $l;

   #$n = shift @a;

    shift @a;
#	shift @a;

    map { $_ = '0.00' if ($_ eq '')} @a;

    #print $n . "\t" . join("\t", @a);
	print join("\t", @a);

    print "\n";
    
}
