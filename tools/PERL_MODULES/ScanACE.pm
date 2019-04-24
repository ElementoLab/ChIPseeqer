package ScanACE;
use lib qw(/home/olly/PERL_MODULES);
use Sets;
use DataFiles;

require "Motif.pm";

my $df = DataFiles->new;

sub new {
    my $self  = {};
    $self->{MOTIF_OBJECT} = undef;
    $self->{SCANACE_MOTIF} = undef;
    $self->{GC} = 0.38;
    $self->{SCANACE} = $df->get("SCANACE");
    $self->{RUN} = 1;
    $self->{THRESHOLD} = undef;
    $self->{RESFILE}   = undef;
    $self->{MOTIF_AVERAGE}   = undef;
    $self->{MOTIF_STDDEV}    = undef;
    $self->{MOTIFS}    = [];

    $self->{VERBOSE} = 0;
    $self->{NBMOTIFS} = undef;
    $self->{STDDEV} = undef;
    $self->{SITES} = [];
    
    bless($self);           # but see below
    return $self;
}

sub setVerbose {
    my ($self, $s)  = @_;
    $self->{VERBOSE} = $n;
}


sub setMotif {
    my ($self, $s_motifFile)  = @_;
    $self->{SCANACE_MOTIF} = $s_motifFile;
}

sub setMotifObject {
    my ($self, $s_motifFile)  = @_;
    
    $self->{MOTIF_OBJECT} = Motif->new;
    $self->{MOTIF_OBJECT}->readFromMotifFile($s_motifFile);
    
    # should output to ScanACE
}


#
#  read file of motifs
#
sub readAlignACEMotifs {
  my ($self, $file) = @_;
  open IN, $file or die "toto\n";
  my $motifdata = 0;
  my $m         = undef;

  while (my $l = <IN>) {
    chomp $l;
    next if ($l eq "");
    if (($motifdata == 0) && ($l =~ /Motif/)) {
      $motifdata = 1;
    }
    if ($motifdata == 1) {
      if ($l =~ /Motif/) {
	if (defined($m)) {
	  push @{ $self->{MOTIFS} }, $m;
	} 
	$m = "Motif 1\n";
      } elsif ($l =~ /MAP/) {
	# nothing ?
      } else {
	my @a = split /\t/, $l, -1;
	$m .= "$a[0]\n";
      }
    }
    
  } 
  close IN;
}

sub getMotifsTXT {
  my ($self) = @_;
  
  return $self->{MOTIFS};
  
}


sub getNbMotifsInFile {
    my ($self, $s_motifFile)  = @_;
    
    open IN, $s_motifFile or die "Impossible to open $s_motifFile\n" ;
    my @a  = <IN>;
    close IN;
    
    my $i =0;
    foreach my $l (@a) {
	#print $l;
	$i++ if ($l =~ /Motif/);
    }

    return $i;
}

sub setGC {
    
    my ($self, $f_gc)  = @_;

    $self->{GC} = $f_gc;
    
}

sub setThresholdScore {
    
    my ($self, $f_thr)  = @_;

    $self->{THRESHOLD} = $f_thr;
    
    
}


sub setNoThresholdScore {
    
    my ($self)  = @_;

    $self->{THRESHOLD} = undef;
    
    
}

sub setNbMotifs {

    my ($self, $n)  = @_;
    
    $self->{NBMOTIFS} = $n;
    
}

sub setFasta {
    my ($self, $n)  = @_;

    $self->{FASTAFILE} = $n;
}


sub setStdDev {
    
     my ($self, $n)  = @_;

     # we can't have both
     # $self->{THRESHOLD} = undef;
    
     $self->{STDDEV} = $n;
    
}

sub run {
    
    my ($self) = @_;
    
    $self->{RESFILE} = Sets::getTempFile("/tmp/tmp.scanace");

    unlink $self->{RESFILE} . "_n1.scn" if (-e $self->{RESFILE} . "_n1.scn");;


    
    if (!defined($self->{NBMOTIFS}) && !defined($self->{STDDEV})) {
	
	die "Please specify a number of motifs\n";
	
    }


    $s_todo = "$self->{SCANACE} -i $self->{SCANACE_MOTIF} -o $self->{RESFILE} -z $self->{FASTAFILE} -g $self->{GC}"; 
    
    
    if ($self->{NBMOTIFS} > 0) {
	$s_todo .= " -s $self->{NBMOTIFS} ";
    }

    if ($self->{STDDEV} > 0) {
	$s_todo .= " -c $self->{STDDEV} ";
    }
    
    print "$s_todo\n"  if ($self->{VERBOSE} == 1);
    
    if ($self->{RUN}) {
	system(($s_todo)); 
    }
        
    
    # store the results in an array
   

		    
}


sub runScanACE {
    
    my ($self, $s_fastaFile)  = @_;

    my $t_m = $self->{SCANACE_MOTIF};

    # outs a temp motif
    if (!$s_motifFile) {

	my $t = $self->{MOTIF_OBJECT}->getScanAceWM;
	
	$t_m = "/tmp/mmm.ace";
	
	open OUT, ">$t_m";
	print OUT $t;
	close OUT;
    }
    
    $self->{RESFILE} = Sets::getTempFile("/tmp/tmp.scanace");

    unlink $s_resFile . "_n1.scn" if (-e $s_resFile . "_n1.scn");;

    if (!defined($self->{NBMOTIFS}) && !defined($self->{STDDEV})) {
	
	die "Please specify a number of motifs\n";
	
    }


    $s_todo = "$self->{SCANACE} -i $t_m -o $self->{RESFILE} -z $s_fastaFile -g $self->{GC}"; 
    if ($self->{NBMOTIFS} > 0) {
	$s_todo .= " -s $self->{NBMOTIFS} ";
    }

    if ($self->{STDDEV} > 0) {
	$s_todo .= " -c $self->{STDDEV} ";
    }
    
    print "$s_todo\n"; # if ($self->{VERBOSE} == 1);
    
    if ($self->{RUN}) {
	system(($s_todo)); 
    }
        
    
    
    unlink "/tmp/mmm.ace";

    # store the results in an array
   
    # $self->_readScanaceResults($self->{RESFILE});
		    
}

sub clean {
    my ($self)  = @_;

    system("rm $self->{RESFILE}" . "_n*.scn");
}

sub setRun {
    my ($self, $r)  = @_;
 
    $self->{RUN} = $r;

}





sub getSites {
    
    my ($self, $n) = @_;

    if (!$n) {
	$n = 1;
    }
    
    $self->{SITES} =  $self->_readScanaceResults($self->{RESFILE} . "_n$n.scn");

    return $self->{SITES};
}



sub getBestSites {
    
    my ($self, $n) = @_;

    if (!$n) {
	$n = 1;
    }
    
    $self->{SITES} =  $self->_readScanaceResults($self->{RESFILE} . "_n$n.scn");

    my %H = ();

    foreach my $r (@{ $self->{SITES} }) {
	if (!defined($H{$r->[0]}) || ($r->[4] > $H{$r->[0]}->[4])) {
	    $H{$r->[0]} = $r;
	}  
    }

    my @aa = values(%H);

    

    @aa = sort { $b->[4] <=> $a->[4] } @aa;

    #foreach my $r (@aa) {
#	print join("\t", @$r) . "\n";
#    }

    return \@aa;
}


sub printResFile {
	my ($self, $n) = @_;

	if (!$n) {
        $n = 1;
	    }

	system("cat $self->{RESFILE}_n$n.scn");
	
}


sub getUniqueORFAboveThreshold {
    
    my ($self, $n) = @_;
    
    if (!$n) {
	$n = 1;
    }
    
    $self->{SITES} =  $self->_readScanaceResults($self->{RESFILE} . "_n$n.scn");
   
         
     my @a_tmp = ();

     my %index = ();

     foreach my $s (@{$self->{SITES}}) {
	 
	 if ($index{$s->[0]} != 1) {
	     push @a_tmp, $s->[0];
	     $index{$s->[0]} = 1;
	 }
	 
     }
     
     return \@a_tmp;
     
} 

sub getAverage {
    my ($self) = @_;
    return  $self->{MOTIF_AVERAGE};
}


sub getStdDev {
    my ($self) = @_;
    return  $self->{MOTIF_STDDEV};
}





sub _readScanaceResults {
    my ($self, $s_file) = @_;
     open IN, $s_file or die "Could not open ScanAce file\n";
     
     my @a_motifs = ();
     my $s_flag = 0;
     while (my $s_line = <IN>) {

	 #print $s_line;
	 chomp $s_line;
	 next if ($s_line eq "");

	 if ($s_line =~ /Motif Average\: ([\d\.]+)\tStd\. Dev\.\: ([\d\.]+)/) {

	     $self->{MOTIF_AVERAGE} = $1;
	     $self->{MOTIF_STDDEV}  = $2;

	     #print "t=$1\t$2\n";
	     #print "\n";
	     #<STDIN>;
	 }


	 if ($s_line =~ /Sites:/) {
	     $s_flag = 1;
	     next;
	 }	 
	 if ($s_flag == 1) {

	     my $pattern = 'upstream';
	     next if ($s_line =~ /$pattern/);
	     
	     
	     $pattern = 'end';
	     next if ($s_line =~ /$pattern/);
	     
	     #print "l=" . $s_line . "\n";
	     
	     my @a_tmp = split /\t/, $s_line; 

	     #print "OK, score= $a_tmp[4]\n";


	     if (defined($self->{THRESHOLD}) && ( $self->{THRESHOLD} > $a_tmp[4])) {
		 #print "stop because line = $s_line\n";

		 
		 last;
	     }
	     
	     #print "OK BIS, score= $a_tmp[4]\n";

	     next if ($a_tmp[4] < 0.0);
	     
	     push @a_motifs, \@a_tmp;
	     
	 } else {
	     next;
	 }
     }     
     close IN;
    #exit;

    #print scalar(@a_motifs) . " entries ..\n";

     return \@a_motifs;
}




1;
