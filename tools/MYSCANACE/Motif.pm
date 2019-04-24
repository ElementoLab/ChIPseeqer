package Motif;
use Sets;
#use strict;

sub new {
    my $self  = {};
    $self->{LENGTH}        = 0;
    $self->{HEIGHT}        = 0;
    $self->{STARS}         = [];
    $self->{MOTIF}         = [];
    $self->{BACKGROUND}    = [];
    $self->{DIBACKGROUND}  = [];
    $self->{VERBOSE}       = 0;
    $self->{NAME}          = undef;
    $self->{MARKED}        = 0;
    $self->{SCORE}      = undef;
    $self->{ SEQNAMES_FOR_SITE } = [];
    $self->{STAR_STR}      = undef;
    srand;
    $self->{SITES}         = []; # contains all the sites
    bless($self);                # but see below
    return $self;
}

#
# return the number of symbols that have more than X symbols
#
sub getNbColumnsWithMaxNbSymbols {
    my ($self, $n) = @_;
    
    my $cnt = 0;
    #print $self->{LENGTH} . "\n";
    for (my $i=0; $i<$self->{LENGTH}; $i++) {
	my $notzeros = 0;
	foreach my $nt (("A", "C", "T", "G")) {
	    #print $nt . "\t" . $self->{MOTIF}->[$i]->{$nt} . "\n";
	    if ($self->{MOTIF}->[$i]->{$nt} > 0.0) {
		$notzeros ++;
	    }
	}
	if ($notzeros <= $n) {
	    $cnt ++;
	} 

    }
    return $cnt;
}

sub mark {
     my ($self, $n) = @_;
    $self->{MARKED}  = $n;
}

sub setScore {
    my ($self, $n) = @_;
    $self->{SCORE}  = $n;
}

sub getScore {
    my ($self) = @_;
    return $self->{SCORE}; 
}


sub setName {
    my ($self, $n) = @_;
    $self->{NAME}  = $n;
}


sub getName {
    my ($self, $n) = @_;
    return $self->{NAME};
}




# recalc the WM from the stored sites
sub reCalc {


}


sub getSites {
    my ($self) = @_;
    
    return $self->{SITES};
}

sub getStars {
    my ($self) = @_;
    
    return $self->{STARS};
}

# return the size of a motif
sub getSize {
    
    my ($self) = @_;
    
    return $self->{LENGTH};
    
}


sub readFromMotifFile {
    
    my ($self, $s_motifFile)  = @_;

    
    open MOT, $s_motifFile;
    my @a_motif = ();
    # read the consensus line
    my $s_line = <MOT>;

    my $s_consensus = $s_line;
    
    chomp $s_line;
    my @a_tmp = split //, $s_line;
    $self->{LENGTH} = length($s_line);
    
    for (my $i=0; $i<scalar(@a_tmp); $i++) {
	if ($a_tmp[$i] eq "*") {
	    $self->{STARS}->[$i] = 1;	
	} else {
	    $self->{STARS}->[$i] = 0;	
	}
    }
    
    # init the motif counts
    for (my $i=0; $i<scalar(@a_tmp); $i++) {
	$self->{MOTIF}->[$i]->{A} = 0;
	$self->{MOTIF}->[$i]->{C} = 0;
	$self->{MOTIF}->[$i]->{T} = 0;
	$self->{MOTIF}->[$i]->{G} = 0;

	$self->{MOTIF}->[$i]->{AA} = 0;
	$self->{MOTIF}->[$i]->{CA} = 0;
	$self->{MOTIF}->[$i]->{TA} = 0;
	$self->{MOTIF}->[$i]->{GA} = 0;
	$self->{MOTIF}->[$i]->{AC} = 0;
	$self->{MOTIF}->[$i]->{CC} = 0;
	$self->{MOTIF}->[$i]->{TC} = 0;
	$self->{MOTIF}->[$i]->{GC} = 0;
	$self->{MOTIF}->[$i]->{AT} = 0;
	$self->{MOTIF}->[$i]->{CT} = 0;
	$self->{MOTIF}->[$i]->{TT} = 0;
	$self->{MOTIF}->[$i]->{GT} = 0;
	$self->{MOTIF}->[$i]->{AG} = 0;
	$self->{MOTIF}->[$i]->{CG} = 0;
	$self->{MOTIF}->[$i]->{TG} = 0;
	$self->{MOTIF}->[$i]->{GG} = 0;
	
    }
    
    $self->{HEIGHT} = 0;

    while ($s_line = <MOT>) {    
	chomp $s_line;
	next if ($s_line eq '');
	my @a_tmp = split //, $s_line;
	for (my $i=0; $i<scalar(@a_tmp); $i++) {
	    $self->{MOTIF}->[$i]->{ $a_tmp[$i] } ++;

	    # dinucl stuff
	    if ($i>0) {
		
		$self->{MOTIF}->[$i]->{ $a_tmp[$i-1] . $a_tmp[$i] } ++;

	    }
	}
	$self->{HEIGHT}++;	

	# also store the site
	push @{$self->{SITES}}, $s_line;
    }    
}



#
# reads a motif from a ScanACE WM (simple version, one motif / file)
#
sub readScanACEMotif {
    
    my ($self, $s_motifFile)  = @_;

    
    open MOT, $s_motifFile or die "Cannot read that matrix ($s_motifFile ..\n";
    
    my @a_motif = ();
    
    # read all the motif
    my @a_lines = <MOT>;
    chomp @a_lines;

    close MOT;

    # remove the first line
    shift @a_lines;
    
    # the second line should contain one site
    my @a       = split /\t/, $a_lines[0];
    $s_line     = $a[0];
    my $i_width = length($s_line);
    
    # the width of the motif is 0 !??
    die "motif width == 0 ..\n" if ($i_width == 0);
    
    for (my $i=0; $i<$i_width; $i++) {
	$self->{MOTIF}->[$i]->{A} = 0;
	$self->{MOTIF}->[$i]->{C} = 0;
	$self->{MOTIF}->[$i]->{T} = 0;
	$self->{MOTIF}->[$i]->{G} = 0;

	$self->{MOTIF}->[$i]->{AA} = 0;
	$self->{MOTIF}->[$i]->{CA} = 0;
	$self->{MOTIF}->[$i]->{TA} = 0;
	$self->{MOTIF}->[$i]->{GA} = 0;
	$self->{MOTIF}->[$i]->{AC} = 0;
	$self->{MOTIF}->[$i]->{CC} = 0;
	$self->{MOTIF}->[$i]->{TC} = 0;
	$self->{MOTIF}->[$i]->{GC} = 0;
	$self->{MOTIF}->[$i]->{AT} = 0;
	$self->{MOTIF}->[$i]->{CT} = 0;
	$self->{MOTIF}->[$i]->{TT} = 0;
	$self->{MOTIF}->[$i]->{GT} = 0;
	$self->{MOTIF}->[$i]->{AG} = 0;
	$self->{MOTIF}->[$i]->{CG} = 0;
	$self->{MOTIF}->[$i]->{TG} = 0;
	$self->{MOTIF}->[$i]->{GG} = 0;
    }
    
    
    $self->{LENGTH} = $i_width;
    

    foreach my $s_line (@a_lines) {
	# empty line
	next if ($s_line eq '');
	
	# kept the left part
	my @a = split /\t/, $s_line;
	$s_line = $a[0];
	

	my @a_tmp = split //, $s_line;
	for (my $i=0; $i<scalar(@a_tmp); $i++) {
	    # ****
	    if ($s_line =~ /\*/) {
		
		if ($a_tmp[$i] eq "*") {
		    $self->{STARS}->[$i] = 1;	
		} else {
		    $self->{STARS}->[$i] = 0;	
		}
		
	    } else {
		$self->{MOTIF}->[$i]->{ $a_tmp[$i] } ++;
		
		# dinucl stuff
		if ($i>0) {
		    $self->{MOTIF}->[$i]->{ $a_tmp[$i-1] . $a_tmp[$i] } ++;
		}
	    }
	    
	}
	 
	if ($s_line !~ /\*/) {
	    $self->{HEIGHT}++;	
	    # also store the site
	    push @{$self->{SITES}}, $s_line;
	}
	
    }    
}

#
# 
#
sub readClassicalWM {
    my ($self, $file) =@_;

    open IN, $file or die "Cannot open $file\n";
    
    my @nt = ("A", "C", "T", "G");
    
    $self->{HEIGHT} = 1;
    my $j = 0;
    while (my $l = <IN>) {

	chomp $l;
		
	
	print "Read $l\n";
	
	
	my @a = split /\t/, $l;
	
	$self->{LENGTH} = scalar(@a);
	
	for (my $i=0; $i<$self->{LENGTH}; $i++) {
	    $self->{MOTIF}->[$i]->{$nt[$j]} = $a[$i];
	}	
	
	$j++;
	
    }
    
    close IN;
    
    for (my $i=0; $i<$self->{LENGTH}; $i++) {
	$self->{STARS}->[$i] = 1;
    }
}



#
# 
#
sub readPollardWM {
    my ($self, $file) =@_;

    open IN, $file or die "Cannot open $file\n";
    
    

    my $j = 0;
    while (my $l = <IN>) {
	chomp $l;
	my @a  = split /[\ \t]+/, $l;
	
	my $nt = shift @a;

	$self->{LENGTH} = scalar(@a);
	
	for (my $i=0; $i<$self->{LENGTH}; $i++) {
	    $self->{MOTIF}->[$i]->{$nt} = $a[$i];
	}	
	
	$j++;
	
    }
    
    close IN;

    my $cnt_max = -1;
    for (my $i=0; $i<$self->{LENGTH}; $i++) {
      my $cnt = 0;
      foreach my $k (keys(%{ $self->{MOTIF}->[$i] })) {
	$cnt += $self->{MOTIF}->[$i]->{$k};
      }

      if ($cnt > $cnt_max) {
	$cnt_max = $cnt;
      }
    }
    

    $self->{HEIGHT} = $cnt_max;

    for (my $i=0; $i<$self->{LENGTH}; $i++) {
	$self->{STARS}->[$i] = 1;
    }
}



sub _initCounts {
    my ($self, $n) =@_;
    
    # init the motif counts
    for (my $i=0; $i<$n; $i++) {
	$self->{MOTIF}->[$i]->{A} = 0;
	$self->{MOTIF}->[$i]->{C} = 0;
	$self->{MOTIF}->[$i]->{T} = 0;
	$self->{MOTIF}->[$i]->{G} = 0;
	
	$self->{MOTIF}->[$i]->{AA} = 0;
	$self->{MOTIF}->[$i]->{CA} = 0;
	$self->{MOTIF}->[$i]->{TA} = 0;
	$self->{MOTIF}->[$i]->{GA} = 0;
	$self->{MOTIF}->[$i]->{AC} = 0;
	$self->{MOTIF}->[$i]->{CC} = 0;
	$self->{MOTIF}->[$i]->{TC} = 0;
	$self->{MOTIF}->[$i]->{GC} = 0;
	$self->{MOTIF}->[$i]->{AT} = 0;
	$self->{MOTIF}->[$i]->{CT} = 0;
	$self->{MOTIF}->[$i]->{TT} = 0;
	$self->{MOTIF}->[$i]->{GT} = 0;
	$self->{MOTIF}->[$i]->{AG} = 0;
	$self->{MOTIF}->[$i]->{CG} = 0;
	$self->{MOTIF}->[$i]->{TG} = 0;
	$self->{MOTIF}->[$i]->{GG} = 0;
	
    }
    
}


#
# create one motif from a single site
#
sub createMotifFromSite {
    
    my ($self, $s, $stars, $n) = @_;
    
    if (!defined($n)) {
	$n = 10;
    }
    
    
    $self->_initCounts(length($s));
    for (my $i=0; $i<$n; $i++) {
	$self->addSite($s);
    }
    $self->setStars($stars);
    
}


sub setStars {
    my ($self, $s) = @_;

    #print "$s\n";
    $self->{STAR_STR} = $s;

    my @a_tmp = split //, $s;
    for (my $i=0; $i<scalar(@a_tmp); $i++) {
	if ($a_tmp[$i] eq "*") {
	    $self->{STARS}->[$i] = 1;	
	} else {
	    $self->{STARS}->[$i] = 0;	
	}

    }
    
}


sub addSiteFromSeq {
    my ($self, $s, $n) = @_;
    $self->addSite($s);
    push @{ $self->{ SEQNAMES_FOR_SITE } }, $n;
}

sub getNumberDistinctSequences {
    my ($self, $n) = @_;
    
    my %H = ();
    foreach my $r (@{ $self->{ SEQNAMES_FOR_SITE } }) {
	$H{ $r } ++;
    }

    return $H{ $n };
}


# add a site $s to a motif
sub addSite {
    
    my ($self, $s) = @_;
    # 
    my @a_tmp = split //, $s;
    for (my $i=0; $i<scalar(@a_tmp); $i++) {
	$self->{MOTIF}->[$i]->{ $a_tmp[$i] } ++;	
	# dinucl stuff
	if ($i>0) {
	    $self->{MOTIF}->[$i]->{ $a_tmp[$i-1] . $a_tmp[$i] } ++;
	}
    }
    $self->{HEIGHT}++;	
    # also store the site
    push @{$self->{SITES}}, $s;

    if (!defined($self->{LENGTH}) || ( $self->{LENGTH} == 0) ) {
	$self->{LENGTH} = length($s);
    }

}


#
#  create a weighted motif
#
sub createMotifFromWeightedSite {
    
    my ($self, $s, $w, $stars) = @_;

    if ($self->{VERBOSE}) {
	print "Creating wm from $s\n";
    }
        
    $self->_initCounts(length($s));    
    $self->addWeightedSite($s, $w);
    $self->setStars($stars);
    $self->{LENGTH} = length($s);

 #   $self->print;
}


#
# add a WEIGHTED site $s to a motif
#
sub addWeightedSite {
    
    my ($self, $s, $w) = @_;
   
    #if ($self->{VERBOSE}) {
    #print " adding $s\n";
    #}
    

    my @a_tmp = split //, $s;

    for (my $i=0; $i<scalar(@a_tmp); $i++) {
	$self->{MOTIF}->[$i]->{ $a_tmp[$i] } += $w;	
    }

    $self->{HEIGHT} += $w;	
      
    #$self->print;
    
}



#
# add a WEIGHTED site $s to a motif
#
sub removeWeightedSite {
    
    my ($self, $s, $w) = @_;
   
    my @a_tmp = split //, $s;

    for (my $i=0; $i<scalar(@a_tmp); $i++) {
	$self->{MOTIF}->[$i]->{ $a_tmp[$i] } -= $w;	
    }

    $self->{HEIGHT} -= $w;	
      
}





# add one random site to the MOTIF
sub addRandomSite {
    
    my ($self) = @_;
   
    my @nt = ("A", "C", "T", "G");
    
    my $s = "";
    for (my $i=0; $i<$self->{LENGTH}; $i++) {
	$s .= $nt[ int(rand(4)) ];
    }

    #print "Adding $s\n";
	
    $self->addSite($s);
    
}


#
# get the ScanACE Wm
#
sub getScanAceWM {

    my ($self) = @_;

    my $s = "Motif 1\n";

    for (my $i=0; $i<$self->{HEIGHT}; $i++) {
	$s .= $self->{SITES}->[$i];
	$s .= "\n";
    }

    for (my $i=0; $i<$self->{LENGTH}; $i++) {
	if ($self->{STARS}->[$i] == 1) {
	    $s .= "*";
	} else {
	    $s .= " ";
	}
    }
    $s .= "\n";
    
    return $s;
    
}

sub printAlignACE {
    my ($self)  = @_;
    
    print $self->getAlignACEString();

}


sub getAlignACEString {
    my ($self)  = @_;

    my $txt = "";
    foreach my $s (@{ $self->{SITES} }) {
	$txt .= "$s\n";
    }
    $txt .= $self->{STAR_STR}; $txt .= "\n";

    return $txt;
}


sub writeSites {
    my ($self, $f)  = @_;
    
    open OUT, ">$f";
    foreach my $s (@{ $self->{SITES} }) {
	print OUT "$s\n";
    }
    close OUT;
}


sub print {
    
    my ($self)  = @_;
    
    print "Here is the motif :\n";
    print "size=$self->{LENGTH}\n";
    print "\t";
    for (my $i=0; $i<$self->{LENGTH}; $i++) {
	print  "*" x $self->{STARS}->[$i] . "\t";
    }
    print "\n";
    foreach my $l (('A', 'C', 'T', 'G')) {
	print "$l\t";
	for (my $i=0; $i<$self->{LENGTH}; $i++) {
	    #print $self->{MOTIF}->[$i]->{$l} . "\t";
	    printf("%5.1f", $self->{MOTIF}->[$i]->{$l});
	}
	print "\n";
    }
    
    
}


sub calcLogo {
    my ($self)  = @_;


    $self->{LOGO} = [];
    
    my $eps     = 0.0001;

    for (my $i=0; $i<$self->{LENGTH}; $i++) {
	my $entropy = 0.0;
	foreach my $l (('A', 'C', 'T', 'G')) {
	    if ($self->{MOTIF}->[$i]->{$l} > 0.0) {
		my $p = $self->{MOTIF}->[$i]->{$l} / $self->{HEIGHT};
		$entropy += $p * log($p) / log(2.0);
	    }
	}

	$entropy = -$entropy;
	$R       = 2.0 - ($entropy + $eps);
	
	foreach my $l (('A', 'C', 'T', 'G')) {
	    my $p = $self->{MOTIF}->[$i]->{$l} / $self->{HEIGHT};
	    $self->{LOGO}->[$i]->{$l} = $p * $R;

	    #print "$i $l --> $p * $R\n"; 
	}
	
    }

    return $self->{LOGO};

}


sub printSimple {
    
    my ($self)  = @_;
  
    #for (my $i=0; $i<$self->{LENGTH}; $i++) {
#	print  $self->{STARS}->[$i] . "\t";
    #}
    #print "\n";
    
    foreach my $l (('A', 'C', 'T', 'G')) {
	
	
	for (my $i=0; $i<$self->{LENGTH}; $i++) {
	    next if ($self->{STARS}->[$i] == 0);
	    #print "$l\t";
	    
	    printf("%3.2f\t", $self->{MOTIF}->[$i]->{$l}/$self->{HEIGHT});
	}
	print "\n";
    }
    
    
}


sub getSimpleMatrix {
    
    my ($self)  = @_;
      
    my $s = "";

    foreach my $l (('A', 'C', 'T', 'G')) {
		
	for (my $i=0; $i<$self->{LENGTH}; $i++) {

	    if ($self->{STARS}->[$i] == 0) {
		$s .= sprintf("%3.2f\t", -1.0);
	    } else {
		$s .= sprintf("%3.2f\t", $self->{MOTIF}->[$i]->{$l}/$self->{HEIGHT});
	    }
	}
	$s .= "\n";
    }
    
    return $s;
    
}



sub writeSimpleMatrix {
    my ($self, $f) = @_;
    
    my $txt = $self->getSimpleMatrix();
    open OUT, ">$f" or die "cannot open file $f\n";
    print OUT $txt;
    close OUT;
}




sub readBackgroundFromFile {
    my ($self, $s_file) = @_;
    
    open DB, $s_file or die "Could not open BKG file\n";
    my %h_bkg = ();
    while (my $s_line = <DB>) {
	chomp $s_line;
	my @a_tmp = split /\t/, $s_line; 
	$h_bkg{ $a_tmp[0] } = $a_tmp[1];
	#print "h{$a_tmp[0]} = $a_tmp[1]\n";
    }
    $self->{BACKGROUND} = \%h_bkg;
}


sub setBkgFrequencies {
    
    my ($self, $h_ref1, $h_ref2) = @_;
    
    $self->{BACKGROUND}   = $h_ref1;
    $self->{DIBACKGROUND} = $h_ref2;
    
}

#
#  set manual frequencies
#
#
sub setManualBkgFrequencies {
    my ($self, $a, $c, $t, $g) = @_;
    
    my %a = ( "A" => $a,
	      "C" => $c,
	      "T" => $t,
	      "G" => $g );

    $self->{BACKGROUND}   = \%a;
    $self->{DIBACKGROUND} = undef;
    
}



sub printBkg {
    my ($self) = @_;

    foreach my $l (('A', 'C', 'T', 'G')) {
	print "BKG{$l} = $self->{BACKGROUND}->{$l}\n";
    }
}


sub _diFreq {
    
    my ($self, $i, $j, $k, $l) = @_;

    my $sum = 0;

    for (my $a=0; $a<$self->{HEIGHT}; $a++) {
	
	my @t = split //, $self->{SITES}->[$a];

	$sum += 1 if (($t[$i] eq $k) && ($t[$j] eq $l));
	
    } 

    #print "sum of joint occurences $k $l = $sum\n";
    
    #$sum += 0.00000001;
    

    return $sum / $self->{HEIGHT};
    
}


# calculate the mutual information between columns of the alignment (sites)
sub mutualInformationBetweenCols {
    my ($self) = @_;

   
    
    for (my $i=0; $i<$self->{LENGTH}-1; $i++) {
	#for (my $j=$i; $j<$self->{LENGTH}; $j++) {
	my $j = $i + 1;
	    my $mut = 0.0;

	    foreach my $k (('A', 'C', 'T', 'G')) {
		
		
		foreach my $l (('A', 'C', 'T', 'G')) {
		
		    # how many times do we have ($k, $l) togetherm $k in seq1, $l in seq2 ? 
		    my $XY = $self->_diFreq($i, $j, $k, $l);
		    
		    #print "Sum of occurences of $k in col $i = " . $self->{MOTIF}->[$i]->{$k} . "\n";
		    my $X  = ($self->{MOTIF}->[$i]->{$k} + 0.00000001) / $self->{HEIGHT};
		    #print "Sum of occurences of $l in col $j = " . $self->{MOTIF}->[$j]->{$l} . "\n";
		    my $Y  = ($self->{MOTIF}->[$j]->{$l} + 0.00000001) / $self->{HEIGHT};
		    
		    next if ($XY == 0);
		  
		    $mut += $XY * $self->_log2 ( $XY / ( $X * $Y ) ); 

		}

	    }

	    print "Mutual information between col $i and $j = $mut\n";

	#}
    }
    
}


sub getKL {
    my ($self) = @_;
    
    my $sum = 0.0;
    
    for (my $j=0; $j<$self->{LENGTH}; $j++) {
	next if ($self->{STARS}->[$j] == 0);
	foreach my $i (('A', 'C', 'T', 'G')) {
	    
	    my $f = ($self->{MOTIF}->[$j]->{$i} + 0.00000001 )/ $self->{HEIGHT};
	    
	    if ($self->{VERBOSE}) {

		 print sprintf("kl [$i] += ( %3d / %3d ) * log2 ( [%3d / %3d] / %3.2f = %f\n",
			  ($self->{MOTIF}->[$j]->{$i} + 0.00000001 ), $self->{HEIGHT},
			       ($self->{MOTIF}->[$j]->{$i} + 0.00000001 ), $self->{HEIGHT},
			       $self->{BACKGROUND}->{$i},
			       $f * $self->_log2( $f / $self->{BACKGROUND}->{$i} ));
		} 
	    $sum += $f * $self->_log2( $f / $self->{BACKGROUND}->{$i} );
	}
    }
    
    if ($self->{VERBOSE}) {
	print "kl = $sum\n";
    }

    return $sum;
}

#
#  compare the current motif with an entire library
#
sub compareACELibrary {


}

sub compareACEScoreWithMotif {
    my ($self, $other) = @_;

    my $m1tmp = "/tmp/m1.ace";
    open TMP, ">$m1tmp";
    print TMP $self->getScanAceWM;
    close TMP;

    my $m2tmp = "/tmp/m2.ace";
    open TMP, ">$m2tmp";
    print TMP $other->getScanAceWM;
    close TMP;
    
    $s_todo = "CompareACE $m1tmp $m2tmp";
    $score = `$s_todo`;
    chomp $score;

    #unlink $m2tmp;
    #unlink $m1tmp;
    
    return $score;
    
}

sub KLdivergenceWithMotif {
    my ($self, $other) = @_;
    
    my $sum = 0.0;
    
    for (my $j=0; $j<$self->{LENGTH}; $j++) {
	next if ($self->{STARS}->[$j] == 0);
	foreach my $i (('A', 'C', 'T', 'G')) {
	    my $p = ($self->{MOTIF}->[$j]->{$i} + 0.00000001 )/ $self->{HEIGHT};
	    my $q = ($other->{MOTIF}->[$j]->{$i} + 0.00000001 )/ $other->{HEIGHT};
	    #print sprintf("f[%d,%s]=%5.4f\n", $j, $i, $f);
	    $sum += $p * log( $p / $q );
	}
    }

    return $sum;
}

sub _log2 {

    my ($self, $x) = @_;

    return log($x) / log(2);
    
}



# calculate the score of a binding site (1rst order M model)
sub score1M {

    my ($self, $s_seqtoscore) = @_;

    #print "get the score of $s_seqtoscore\n";

    my $N = $self->{HEIGHT};
    my @a_seq = split //, $s_seqtoscore;

    my $S = 0.0;

    # initialize with the first 1nt 
    my $b = $a_seq[0];
    
    #print $self->{MOTIF}->[0]->{$b};
    #print $self->{BACKGROUND};
    #print "\n$b\n";
    #<STDIN>;

    $S   += $self->_log2 ( ( $self->{MOTIF}->[0]->{$b} + 
		    $self->{BACKGROUND}->{$b} ) / ( $N + 1 ) ) - log ( $self->{BACKGROUND}->{$b} );


    # now goes thru the whole seq
    for (my $i=1; $i<$self->{LENGTH}; $i++) {
	
	#get the current dint (i-1, i) 
        my $a = $a_seq[$i-1]; 
	my $b = $a_seq[$i];

	#print "S += log ( ( self->{MOTIF}->[$i]->{ $a . $b } + self->{DIBACKGROUND}->{$a$b}) / ( $N + 1 ) ) - log ( self->{DIBACKGROUND}->{$a$b} )\n";

	
        $S += $self->_log2 ( ( $self->{MOTIF}->[$i]->{ $a . $b } + $self->{DIBACKGROUND}->{$a . $b}) / ( $N + 1 ) ) - log ( $self->{DIBACKGROUND}->{$a . $b} );
    }

    return $S;
}




# calculate the CLASSICAL score of a binding site (1rst order M model)
sub score0M {

    my ($self, $s_seqtoscore) = @_;

    #print "get the score of $s_seqtoscore\n";

    my $N = $self->{HEIGHT};
    my @a_seq = split //, $s_seqtoscore;

    my $S = 0.0;


    # now goes thru the whole seq
    for (my $i=0; $i<$self->{LENGTH}; $i++) {
	
	my $b = $a_seq[$i];
	
	$S   += $self->_log2 ( ( $self->{MOTIF}->[$i]->{$b} + $self->{BACKGROUND}->{$b} ) / ( $N + 1 ) ) - log ( $self->{BACKGROUND}->{$b} );

    }

    return $S;
}



# calculate the CLASSICAL score of a binding site (1rst order M model)
sub score0M_simple {

    my ($self, $s_seqtoscore) = @_;

    #print "get the score of $s_seqtoscore\n";

    my $N = $self->{HEIGHT};
    my @a_seq = split //, $s_seqtoscore;

    my $S = 1.0;


    # now goes thru the whole seq
    for (my $i=0; $i<$self->{LENGTH}; $i++) {
	
	my $b = $a_seq[$i];
	
	$S   *=  ( $self->{MOTIF}->[$i]->{$b} + $self->{BACKGROUND}->{$b} ) / ( $N + 1 ) ;

    }

    #print "$s_seqtoscore\t$S\n";

    return $S;
}


sub getThresholdByStdDev {
    my ($self, $n) = @_;

    my $avg = 0;
    my $std = 0;
    my @sco = ();

    # for each seq, calc a score
    my $cnt = 0;
    foreach my $s (@{$self->{SITES}}) {
	my $s1 = $self->score0M_simple($s);
	push @sco, $s1;
	$avg += $s1;
	$cnt++;
	
	#my $s2 = $self->score0M_simple(Sets::getComplement($s));
	#push @sco, $s2;
	#$avg += $s1;
	#$cnt++;
    }

    
    $avg /= $cnt;
    
    
    for (my $i=0; $i<$cnt; $i++) {
	$std += ($sco[$i] - $avg) * ($sco[$i] - $avg);
    }
    $std /= ($cnt - 1);
    $std = sqrt($std);
 
    return $avg - $n * $std;
}




# returns an (sorted) array of (position, score, strand, sequence)
#    up to the given threshold
sub scanSequence {
    my ($self, $seq, $t) = @_;

    my @a_seqs = ();

    for (my $i=0; $i<length($seq) - $self->{LENGTH}; $i++ ) {

	my $ss = substr($seq, $i, $self->{LENGTH});
	my $sc = Sets::getComplement($ss);
	my $score1 = $self->score0M_simple($ss);
	my $score2 = $self->score0M_simple($sc);
	
	if ($score1 > $t) {
	    my @a_tmp = ($i, $score1, 0, $ss);
	    push @a_seqs, \@a_tmp;
	}

	if ($score2 > $t) {
	    my @a_tmp = ($i, $score2, 1, $sc);
	    push @a_seqs, \@a_tmp;
	}

    }
    
    
    my @a_sorted = sort { $b->[1] <=> $a->[1] } @a_seqs;
    
    return \@a_sorted;
}


sub generateRandomSequence {
    my ($self) = @_;

    
    my $s = "";
    
    for (my $j=0; $j<$self->{LENGTH}; $j++) {
	
	my @a_cum = ();
	my $i     = 1;
	my @a     = ('A', 'C', 'T', 'G');
	
	$a_cum[0] = 0.0;
	foreach my $l1 (@a) {
	    $a_cum[$i] = $a_cum[$i - 1] + $self->{MOTIF}->[$j]->{$l1} / $self->{HEIGHT};
	    $i++;
	}

	$a_cum[$i - 1] = 1.0;
    
	my $d = rand;

	my $n = undef;
	
	for ($i=1; $i<=4; $i++) {
	    if (($d > $a_cum[$i-1]) && ($d <= $a_cum[$i])) {
		$n = $a[$i-1];
		last;
	    }
	}
	
	
	$s .= $n;
    }

    return $s;

}


sub generateScanACEMotif {
    my ($self)  = @_;
    
    my $s = "Motif 1\n";
    for (my $i=0; $i<100; $i++) {
	$s .= $self->generateRandomSequence;
	$s .= "\n";
    }
    $s .= "*" x $self->{LENGTH};
    $s .= "\n";
    return $s;

}


sub normalize {
    my ($self)  = @_;
    # get a list of position in an array
    
    my @a = ();
    for (my $i=0; $i<$self->{LENGTH}; $i++) {
	$self->{MOTIF}->[$i]->{A} /= $self->{HEIGHT};
	$self->{MOTIF}->[$i]->{C} /= $self->{HEIGHT};
	$self->{MOTIF}->[$i]->{T} /= $self->{HEIGHT};
	$self->{MOTIF}->[$i]->{G} /= $self->{HEIGHT};
    }

}

#
#  returns whether a motif is a palindrome 
#
sub getPalindromeScore {
    my ($self, $wmin, $wmax, $ks, $dirs)  = @_;
    # get a list of position in an array

    #$self->normalize;

    my @a = ();
    for (my $i=0; $i<$self->{LENGTH}; $i++) {
	my @a_tmp =  ($self->{MOTIF}->[$i]->{A},
		      $self->{MOTIF}->[$i]->{C},
		      $self->{MOTIF}->[$i]->{T},
		      $self->{MOTIF}->[$i]->{G});
	push @a, \@a_tmp;
    }

    my $best_k = 0;
    my $best_i = 0;
    my $best_j = 0;
    my $best_c = -10.0;
    my $best_s = "";


    for (my $k=$wmin; $k<=$wmax; $k++) {
	
	for (my $i=0; $i<$self->{LENGTH}-2*$k; $i++) {
	    	    
	    for (my $j=$i+$k; $j<$self->{LENGTH}-$k+1; $j++) {
	    
		# 1. inverted palindromes
		
		# create a first vector
		my @v1 = ();
		my $s1 = 0;
		for (my $l=$i; $l<$i+$k; $l++) {
		    my @a_tmp =  ($self->{MOTIF}->[$l]->{A},
				  $self->{MOTIF}->[$l]->{C},
				  $self->{MOTIF}->[$l]->{T},
				  $self->{MOTIF}->[$l]->{G});
		    push @v1, @a_tmp;
		    $s1++ if ($self->{STARS}->[$l] == 1);
		}


		# create a first vector
		my @v2 = ();
		my $s2 = 0;
		my @v3 = ();
		for (my $l1=$j+$k-1, my $l2=$j; $l1>=$j; $l1--, $l2++) {
		    my @a_tmp =  ($self->{MOTIF}->[$l1]->{T},
				  $self->{MOTIF}->[$l1]->{G},
				  $self->{MOTIF}->[$l1]->{A},
				  $self->{MOTIF}->[$l1]->{C});
		    push @v2, @a_tmp;


		    my @a_tmp =  ($self->{MOTIF}->[$l2]->{A},
				  $self->{MOTIF}->[$l2]->{C},
				  $self->{MOTIF}->[$l2]->{T},
				  $self->{MOTIF}->[$l2]->{G});
		    push @v3, @a_tmp;

		    $s2++ if ($self->{STARS}->[$l2] == 1);

		}
		
		# next if number of signif columns is lesser than window size minus  
		next if (($s1 < $ks) || ($s2 < $ks));

		
		#print "i=$i\tj=$j\n";
		#print join("\t", @v1); print "\n";
		#print join("\t", @v2); print "\n";
		#print join("\t", @v3); print "\n";

		
		#print "$s1\t$s2\t";


		#print Sets::euclidean(\@v1, \@v2);
		
		#print "\t";

		# -1
		if (($dirs == 0) || ($dirs == 2)) {
		    
		    my $c = Sets::pearson(\@v1, \@v2);
		    #print "c(-1)=$c\n";
		    
		    if ($c > $best_c) {
			$best_i = $i;
			$best_j = $j;
			$best_k = $k;
			$best_c = $c;
			$best_s = -1;
		    }
		}
		    
		    
		# 1
		if (($dirs == 1) || ($dirs == 2)) {
		    my $c = Sets::pearson(\@v1, \@v3);
		    #print "c(1)=$c\n";
		    if ($c > $best_c) {
			$best_i = $i;
			$best_j = $j;
			$best_k = $k;
			$best_c = $c;
			$best_s = 1;
		    }
		}
	    
	    }
	    
	}
	
    }

    return  ($best_i, $best_j, $best_k, $best_s, $best_c);
    

}



1;
