# take a snapshot of the files in current dir ..
my $s = `/bin/ls`;
my @a = split /\n/, $s;

my %H = ();
# index
foreach my $l (@a) {
    $H{ $l } = 1;
}


# now indefinite loop
my $n = 1;
while ($n == 1) {

    my $s = `/bin/ls`;
    my @a = split /\n/, $s;
    
    foreach my $r (@a) {
	#print "$r\n";
	if (!defined($H{$r})) {
	    system("echo \"test\" | mail -s \"$r\" elemento\@princeton.edu");
	    print "Email sent for $r\n";
	    $H{$r} = 1;
	} 
    }

    sleep(10);
}
