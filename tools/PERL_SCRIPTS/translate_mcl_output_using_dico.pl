use lib qw(/home/elemento/PERL_MODULES);
use Sets;
use Table;



my $ta = Table->new;
$ta->loadFile($ARGV[1]);
my $h_ref_index = $ta->getIndexKV(0,1);


$ta->loadFile($ARGV[2]);
my $h_ref_desc = $ta->getIndexKV(0,1);


open IN, $ARGV[0];
my @C = ();
my $cnt = 0;
while (my $l = <IN>) {
    chomp $l;
    
    if ($l =~ /^\d/) {
	#print "$l\n";
	my @a = split /[\s]+/, $l;
	shift @a;
	my $n = pop @a;
	push @C, @a;
	
	if ($n eq '$') {
	    my $a_ref_C = Sets::translateSet(\@C, $h_ref_index);
	    &show_cluster($a_ref_C, $cnt++, $h_ref_desc);
	    @C = ();
	} else {
	    push @C, $n;
	}
	
    


    } elsif ($l =~ /^ /) {
	#print "$l\n";
		
	my @a = split /[\s]+/, $l;
	my $n = pop @a;
	shift @a;
	push @C, @a;

	if ($n eq '$') {
	    my $a_ref_C = Sets::translateSet(\@C, $h_ref_index);
	    &show_cluster($a_ref_C, $cnt++, $h_ref_desc);
	    @C = ();
	} else {
	    push @C, $n;
	}
	
    }
    
}

close IN;


sub show_cluster {
    
    my ($r, $n, $desc) = @_;
    
    my $m = scalar(@$r);
    print "CLUSTER $n has $m members\n";
    #print join("\n", @$r); print "\n";

    foreach my $o (@$r) {
	print $o;

	if (defined($desc)) {
	    print "\t". $desc->{$o};
	}

	print "\n";
    }

}
