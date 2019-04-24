package Vista;

sub new {
    my ($self) = {};
    $self->{ALIGN_FILE} = undef;
    $self->{SSEQ1} = undef;
    $self->{SSEQ2} = undef;
    $self->{SNUM1} = undef;
    $self->{SNUM2} = undef;
    $self->{SBARS} = undef;
    my %idx = ();
    $self->{KMERS} = \%idx;
    
    bless $self;
    return $self;
}

sub readKmers {
    my ($self, $f) = @_;
    open INICH, $f;
    while (my $l = <INICH>) {
	my @a = split /\t/, $l;
	$self->{KMERS}->{$a[0]} = 1;
    }
    close INICH;
    

    for (my $i=0; $i<length($self->{SSEQ1}); $i++) {
	my $kmer = substr($self->{SSEQ1}, $i, 7);
	print "$kmer\n";
	if ($self->{KMERS}->{uc($kmer)} == 1) {
	    substr($self->{SSEQ1}, $i, 7) = uc($kmer);
	}
    }

    print $self->{SSEQ1};
    print "\n";
    
}

sub printAlignment {
     my ($self) = @_;
    for (my $i=0; $i<length($self->{SSEQ1}); $i += 60) {
	
	print substr($self->{SSEQ1}, $i, 60); print "\n";
	print substr($self->{SBARS}, $i, 60); print "\n";
	print substr($self->{SSEQ2}, $i, 60); print "\n";
	print "\n";

    }
}


sub readAlignmentFile {
    my ($self, $f) = @_;
    $self->{ALIGN_FILE} = $f;

    open INICH,  $self->{ALIGN_FILE} or die "cannot open $f\n";

    my @lines = <INICH>;
    chomp @lines;

    close INICH;
    
    my $sseq1 = "";
    my $sseq2 = "";
    my $snum1 = "";
    my $snum2 = "";
    my $sbars = "";

    # 4 lines to recreate
    my $i = 0;
    foreach my $l (@lines) {
	if ($l =~ /^seq1/) {
	    my ($seq) = $l =~ /.{9}(.+)$/;
	    $sseq1 .= $seq;
	} 

	if ($l =~ /^seq2/) {
	    my ($seq) = $l =~ /.{9}(.+)$/;	
	    $sseq2 .= $seq;
	} 

	if ($lines[$i+1] =~ /^seq1/) {
	    if ($l =~ /\s{9}(.+)$/) {
		my $tmp  = $1 . substr(" " x 60, 0, 60 - (length($1)));
		#print "ADDING $tmp*\n";		
		$snum1 .= $tmp;
	    } else {
		#print "ADDING " . (" " x 60) . "*\n";
		$snum1 .= " " x 60;
	    }	    
	    
	}

	if ($lines[$i-1] =~ /^seq2/) {
	    if ($l =~ /\s{9}(.+)$/) {
		my $tmp  = $1 . substr(" " x 60, 0, 60 - (length($1)));

		$snum2 .= $tmp; 
	    } else {
		$snum2 .= " " x 60;
	    }
	}

	if ( ($lines[$i-1] =~ /^seq1/) && ($lines[$i+1] =~ /^seq2/) ) {
	    
	    
	    $sbars .= substr($l, 9);
	}
	
	$i ++;
    }
    
    print "*$snum1*\n";
    print "*$sseq1*\n";
    print "*$sbars*\n";
    print "*$sseq2*\n";
    print "*$snum2*\n";


    $self->{SSEQ1} = $sseq1;
    $self->{SSEQ2} = $sseq2;
    $self->{SNUM1} = $snum1;
    $self->{SNUM2} = $snum2;
    $self->{SBARS} = $sbars;
    

}
1;
