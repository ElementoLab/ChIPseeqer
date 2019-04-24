my @matrix   = ();
my @rownames = ();
my @colnames = ();

loadMatrix($ARGV[0], \@matrix, \@rownames,\@colnames);

my $n     = scalar(@rownames);
my $m     = scalar(@colnames);

my @inverted_matrix = ();
my @a = ();
my @b = ();
for (my $j=0; $j<$m; $j++) {
    if (!defined($ARGV[3])) {
	if ($j < $ARGV[1]) {
	    push @b, $colnames[$j];
	}
	if (($j >= $ARGV[1]) && ($j <$ARGV[2])) {
	    push @a, $colnames[$j];
	}
    } else {
	if ($j < $ARGV[1]) {
	    push @b, $colnames[$j];
	}
	if (($j >= $ARGV[1]) && ($j <$ARGV[2])) {
	    push @a, $colnames[$j];
	}
    }
}
my @inverted_colnames = (@a, @b);
print "\t";
print join("\t", @inverted_colnames); 
print "\n";



for (my $i=0; $i<$n; $i++) {
    
    my @a = ();
    my @b = ();

    for (my $j=0; $j<$m; $j++) {
	if ($j < $ARGV[1]) {
	    push @b, $matrix[$i][$j];
	}
	if (($j >= $ARGV[1])  && ($j <$ARGV[2])) {
	    push @a, $matrix[$i][$j];
	}
    }
    my @a_tmp = (@a, @b);
    print "$rownames[$i]\t";
    print join("\t", @a_tmp); 
    print "\n";
}


sub loadMatrix {
    
    my ($f, $a_ref_matrix, $a_ref_rownames, $a_ref_colnames) = @_;

    open IN, $f;
    my $l = <IN>;
    chomp $l;
    my @a = split /\t/, $l, -1;
    shift @a;
    @{ $a_ref_colnames } = @a;

    while (my $l = <IN>) {
        chomp $l;
        my @a = split /\t/, $l, -1;
        my $r = shift @a;
        push @{ $a_ref_rownames }, $r;
        push @{ $a_ref_matrix }, \@a;
    }
    close IN;
}
