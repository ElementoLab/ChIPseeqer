#
#  transform a Latex tabular table into a tsv output
#


# & => \t
# \tabularnewline => \n
# \\ => \n

while ($l = <STDIN>) {
    
    chomp $l;
    
    print replace($l);
    
}


sub replace {
    my ($l) = @_;

    my $lc = $l;

    my @a_lib = ( [ '&', "\t" ],
		  [ '\\\\tabularnewline', "\n" ],
		  #[ '\\\\\\\\', "\n"],
		  [ '\\\\hline\\ *', ''],
		  [ "\\\\cite\\{.+?\\}", ''],
		  [ '\\$.+?\\$', '']
		  );
    #print "l=$lc ";

    foreach my $r (@a_lib) {
	if ($lc =~ /$r->[0]/) {
	    #print "s/$r->[0]/$r->[1]/g\n";

	    $lc =~ s/$r->[0]/$r->[1]/g;
	}
	#print "next\n";
    }
    #print "becomes $lc\n";
    return $lc;
}
