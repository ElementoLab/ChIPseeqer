#
#  parses a lalign report
#
package lalign;

use lib qw(/home/olly/PERL_MODULES);
use MyBlast;
use strict;
use DataFiles;

my $df = DataFiles->new;

sub new {
    my ($self)  = {};
    $self->{FILE} = undef;
    $self->{MAXEVALUE} = undef;
    $self->{NBALIGNS} = 100;
    $self->{RESULTS} = [];
    bless($self);
    return $self;

}


sub setFile {
    my ($self, $s) = @_;
    $self->{FILE} = $s;
}

sub getResults {
    my ($self, $s) = @_;
    return $self->{RESULTS};
}

sub setMaxEvalue {
     my ($self, $s) = @_;
     $self->{MAXEVALUE} = $s;
}


sub processUsingBlast {
    my ($self, $seq1file, $seq2file) = @_;
    
    my $mb = MyBlast->new;

    #$mb->setVerbose(1);

    $mb->setBlastProgram("blastn");
    
    $mb->setQueryDatabase($seq1file);
    $mb->setDatabaseFile($seq2file);
    $mb->setNbProcessors(2);
    $mb->setEvalueThreshold($self->{MAXEVALUE});
    $mb->setMismatchWeight(-6);
    $mb->setMatchWeight(5);
    $mb->setGapOpening(20);
    $mb->setGapExtension(4);
    $mb->setWordLength(7);
    $mb->setQueryStrand(1);
    
    
    

    #blastall -i dmel.seq -d dort.seq.agam -p blastn -G 30 -E 4 -e 0.0001 -r 5 -q -6 -W 7

    #$mb->setVerbose(1);
    my $a_ref = $mb->blastallUnique;
    
    #<STDIN>;
    
    foreach my $r (@$a_ref) {
	
	my $s1 = $r->{QSEQ};
	my $s2 = $r->{DSEQ};

	my @a1 = split //, $s1;
	my @a2 = split //, $s2;
	
	my $cnt = -1;
	my $mystart = $r->{QFROM};
	my @pos = ();
	for (my $i=0; $i<length($s1); $i++) {
	    $cnt ++ if ($a1[$i] ne '-');
	    if ($a1[$i] eq $a2[$i]) {
		push @pos, $mystart+$cnt;
	    }
	    
	}
	
	
	my %hash = (CONSERVED_POS => \@pos, DSTART => $mystart, DEND => $mystart + $cnt, HSEQ => $s2);
	push @{ $self->{RESULTS} }, \%hash;
    }
    
}


sub setNumberAlignments {
     my ($self, $s) = @_;
     $self->{NBALIGNS} = $s;
}

sub processUsingLalign {
    
    my ($self, $seq1file, $seq2file) = @_;

    my $outfile = Sets::getTempFile("tata");
    my $todo = $df->get("LALIGN") . " -n -f -30 $seq1file $seq2file $self->{NBALIGNS} > $outfile 2> /dev/null";
    #print "$todo\n";
    system($todo);

    #system("cat $outfile");
    
    $self->setFile($outfile);
    $self->{OUTLALIGN} = $outfile;
    $self->process();
}

sub dispose {
    my ($self) = @_;
    
    unlink $self->{OUTLALIGN};
}

sub getLalignOutputFile {
    my ($self) = @_;
    return $self->{OUTLALIGN};
}


sub process {
    my ($self) = @_;
    open IN, $self->{FILE};
    my $l = <IN>; 
    $l = <IN>;
    $l = <IN>;
    $l = <IN>;
    $l = <IN>;

    while ($l = <IN>) {
	if ($l =~ /identity/) {
	    my ($id, $ov, $sc, $ev) = $l =~ /([\d\.]+)\% identity in (\d+) nt overlap\; score\:\ +(\d+) E\(10\,000\)\:\ +([\d\.\+\-e]+)/;

	    #print "$sc\t$ov\t$ev\n";

	    #print $l if !$sc;

	    my $mys1    = "";
	    my $mys2    = "";
	    my $mystart = 0;

	    $l = <IN>;
	    while (1) {
		
		
		$l = <IN>; # numbers
		#print $l;
		
		if ($l =~ /^\-/) {

		    
		    #print "GOT (query seq starts at $mystart) : \n$mys1\n$mys2\n";
		    
		    if ($ev < $self->{MAXEVALUE}) {
			
			my @a1 = split //, $mys1;
			my @a2 = split //, $mys2;
			
			my $cnt = -1;
			
			
			my @pos = ();
			for (my $i=0; $i<length($mys1); $i++) {
			    $cnt ++ if ($a1[$i] ne '-');
			    if ($a1[$i] eq $a2[$i]) {
			    #print "$mystart+$cnt is conserved \n";
				push @pos, $mystart+$cnt;
			    }
			    
			}
			
			
			my $hseq = $mys2; #$hseq =~ s/\-//g;
			my %hash = (CONSERVED_POS => \@pos, DSTART => $mystart, DEND => $mystart + $cnt, HSEQ => $hseq);
			
			push @{ $self->{RESULTS} }, \%hash;
		    
		    }  else {
			#print "$ev >= $self->{MAXEVALUE}\n";
		    }
		    
		    
		    $mys1    = "";
		    $mys2    = "";
		    $mystart = 0;
		    $l = <IN>; last;
		    
		} 
		

		# get the number
		my $lc = $l; $lc =~ s/^\ +//g; my @a = split /\ +/, $lc;
		my $n  = $a[0];
		
		# get the position of the last digit
		$l =~ s/^.{7}//g;
		my @a = split //, $l;

		my $cnt = 0;
		for (my $i=0; $i<length($l); $i++) {
		    if ( ($a[$i] =~ /\d/) && ($a[$i+1] =~ /\ /) ) {
			last;
		    } else {
			$cnt ++;
		    }
		}
		
		my $sta = $n - $cnt;
		#print "$n, $cnt, so seq query match starts at $sta\n";

		
		if ($mys1 eq "") {
		    $mystart = $sta;
		}
		

		$l = <IN>; # seq 
		chomp $l;
		$l =~ s/^.{7}//g;  $mys1 .= $l;
		

		$l = <IN>; # semi-cols

		$l = <IN>; # seq
		chomp $l;
		$l =~ s/^.{7}//g;  $mys2 .= $l;

		$l = <IN>; # numbers
		$l = <IN>; # space

		#$l = <IN>; # end ?
		
		
		#print "ITER FINISHED, mys1 = $mys1\n";
		
	    }

	    
	}
	
	
    }

    

    close IN;

}


1;
