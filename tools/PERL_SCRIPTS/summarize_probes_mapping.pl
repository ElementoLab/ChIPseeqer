open IN, $ARGV[0];
while (my $l = <IN>) {
    chomp $l;    
    my @a = split /\t/, $l;
    $H{ $a[0] }->{ $a[1] } ++ if ($a[3] == 25);
}
close IN;

my @ps = keys(%H);

my %PROBE_SETS = ();
open LOG, ">log_affy_mapping.txt";

foreach my $p (@ps) {

    my @ts = keys(%{ $H{$p} });
    my @ma = ();
    #
    #  for a given probe set, gather all CGs, with the number of probes matching the CG
    #
    foreach my $t (@ts) {
	my $n =  $H{ $p }->{ $t };
	my @tm = ($t, $n); 
	push @ma, \@tm;
    }
    
    #
    #  sort based on the number of probes
    #
    @ma = sort { $b->[1] <=> $a->[1] } @ma;

    my $nbt = scalar(@ma);


    #
    #  assign valid probe sets to genes, only one probe set per gene 
    #
    my %TAKEN = ();
    
    my @txt = ();
    for (my $i=0; $i<$nbt; $i++) {
	
	my $t = $ma[$i][0]; $t =~ s/^.+\_//; # good for -RA .. $t =~ s/\-.+$//g;
	my $n = $ma[$i][1];
	
	
	if (!defined($TAKEN{$t}) && ($n >= 8)) {

	    push @txt, $t;
	    #my @tm = ($p, $n); 
	    #push @{ $PROBE_SETS{ $t } }, \@tm;
	    $TAKEN{ $t } = 1;
	}
    }
    
    if (scalar(@txt) > 0) {
	print "$p\t"; print join("/", @txt); print "\n";
    } else {
	print LOG "$p could not be assigned to any gene ..\n";
    }

}
close LOG;

exit;

my %GENES = ();

foreach my $g (keys(%PROBE_SETS)) {
    my @probes = sort { $b->[1] <=> $a->[1] } @{ $PROBE_SETS{$g} };

    my $p = shift @probes; 
    
    my @tm = ($g, $p->[1]); 
    
    push @{ $GENES{ $p->[0] } }, \@tm;

    #print "$p->[0]\t$g\n";

    foreach my $op (@probes) {
	print LOG "$op->[0]\t$op->[1]\tgene already asssigned\t$p->[0]\t$p->[1]\t$g\n";
    }
    
}

foreach my $p (keys(%GENES)) {
    my @genes = sort { $b->[1] <=> $a->[1] } @{ $GENES{$p} };
    
    my $g = shift @genes; 
    
    print "$p\t$g->[0]\n";
    
    foreach my $og (@genes) {
	print LOG "$og->[0]\t$og->[1]\tprobe already assigned\t$g->[0]\t$g->[1]\t$p\n";
    }
    
    
}

close LOG;
