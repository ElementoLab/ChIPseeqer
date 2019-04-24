package SVMlight;
use strict;

sub new {
    
    my $self  = {};


    $self->{MODEL_FILE}  = "/tmp/svmlight.mdl";
    
    $self->{NBPOS}  = undef;
    $self->{NBNEG}  = undef;

    $self->{CROSSVAL} = 0;

    # array containing several examples
    $self->{EXAMPLES}     = [];
    $self->{RAWEXAMPLES}  = [];
    $self->{WEIGHTS}      = [];
    $self->{AFFECTATIONS} = [];
    $self->{CROSSSVALSIZE} = 1;
    # array containing the class of the examples (-1 or +1)
    $self->{CLASSES}     = [];

    $self->{PREDICTIONS} = [];
    $self->{WRITTEN}     = [];
    $self->{MASK} = [];

    $self->{VERBOSE} = 0;

    $self->{G}    = undef;
    $self->{C}    = undef;
    $self->{K}    = undef;
    $self->{LABELS} = undef;

    $self->{CLASSIFY_ERROR} = "";
    $self->{CLASSIFY_OUTPUT} = "";
    $self->{CLASSIFY_RECALL} = undef;
    $self->{CLASSIFY_PRECISION} = undef;

    bless($self);           # but see below
    return $self;
    
}


sub setCrossvalSize {
    my ($self, $nc) = @_;

    $self->{CROSSSVALSIZE} = $nc;

    srand;

    # does the affectation
    my $n = scalar(@{$self->{RAWEXAMPLES}});
    for (my $i=0; $i<$n; $i++) {
	my $r = int(rand( $nc )) + 1;
	$self->{AFFECTATIONS}->[$i] = $r;

	#if ($self->{VERBOSE} == 1) {
	#    print "affect indiv $i to cluster $r\n";
	#}
    }
    
		   
		   
}



sub getNbPos {
    
    my ($self)  = @_;
    
    return $self->{NBPOS};
    
}


sub getPredictions {
    
    my ($self)  = @_;
    
    return $self->{PREDICTIONS};
    
}




sub getNbNeg {
    
    my ($self)  = @_;
    
    return $self->{NBNEG};
    
}

#
#  learn X n
#
sub learnCrossval {

    my ($self)  = @_;
    
    my $ASE = 0;
    my $ASP = 0;
    my $AER = 0;
    
    # write all examples to disk
    $self->_writeExamples();
    
    for (my $j=1; $j<=$self->{CROSSSVALSIZE}; $j++) {
	
	#
	# create a new file not including file $j
	#
	my $newfile = "/tmp/svmlight.allexa.$j";
	system(">$newfile");
	for (my $k=1; $k<=$self->{CROSSSVALSIZE}; $k++) {
	    if ($k != $j) {
		system("cat /tmp/svmlight.exa.$k >> $newfile");
	    }
	}
	
	#
	#  train the SVM
	#
	$self->_run($newfile, "/tmp/svmlight.mdl.$j");
	
	
	
	#
	#  get labels and classes for this set of examples
	#
	my $n = scalar(@{$self->{RAWEXAMPLES}});
	my $cnt = 0;
	my @a = ();
	my @b = ();
		
	for (my $i=0; $i<$n; $i++) {
	    if ($self->{AFFECTATIONS}->[$i] == $j) {
		for (my $k=0; $k<$self->{WEIGHTS}->[$i]; $k++) {
		    push @a, $self->{CLASSES}->[$i];
		    push @b, $self->{LABELS} ->[$i];
		}
	    }
	}

	#
	#  classify the remaining file
	#
	my $a_ref_pred = $self->_classify("/tmp/svmlight.exa.$j", "/tmp/svmlight.mdl.$j");
	

	#
	#  analyze the predictions  
	#
	my $TP = 0;
	my $TN = 0;
	my $FP = 0;
	my $FN = 0;
	my $n = scalar(@$a_ref_pred);
	for (my $k=0; $k<$n; $k++) {
	    
	    #print $self->{WRITTEN}->[$j]->[$k]->[0] . "\t" . $self->{WRITTEN}->[$j]->[$k]->[1] . "\t" . $a_ref_pred->[$k] . "\n"; 

	    
	    my $p = $self->{WRITTEN}->[$j]->[$k]->[0];
	    my $t = $a_ref_pred->[$k]; 

	    #print "$p\t$t\n";
	    
	    if (($p > 0) && ($t > 0)) {
		$TP ++;
	    } 

	    if (($p > 0) && ($t < 0)) {
		$FP ++; 
	    } 

	    if (($p < 0) && ($t > 0)) {
		$FN ++;
	    } 

	    if (($p < 0) && ($t < 0)) {
		$TN ++;
	    } 
	    
	    
	    
	}
	

	

	my $SE = ($TP==0?0:$TP / ($TP + $FN));
	my $SP = ($TN==0?0:$TN / ($TN + $FP));
	
	my $ER = ($TP + $TN) / ($TP + $TN + $FN + $FP);
	

	print "TP = $TP, FN = $FN, TN = $TN, FP = $FP, "; print "SE = $SE, SP = $SP, ER = $ER\n";

	$ASE += $SE;
	$ASP += $SP;
	$AER += $ER;
	
    }

    $ASE /= $self->{CROSSSVALSIZE};
    $ASP /= $self->{CROSSSVALSIZE};
    $AER /= $self->{CROSSSVALSIZE};
    
    print "ASE = $ASE, ASP = $ASP, AER = $AER\n";

}


sub _run {
    my ($self, $exa, $mdl)  = @_;

    my $home = `echo \$HOME`; chomp $home;


    my $s_todo  = "$home/PERL_MODULES/PROGRAMS/SVMLIGHT/svm_learn ";
	
    $s_todo .= " -t 0 " if (!defined($self->{K}) || ($self->{K} == 0));
    $s_todo .= " -t $self->{K} " if ($self->{K});
    $s_todo .= " -g $self->{G} " if ($self->{G});
    $s_todo .= " -c $self->{C} " if ($self->{C});
    
	#$s_todo .= " -g $self->{G} " if (defined $self->{G});
    $s_todo .= " -x 1 " if ($self->{CROSSVAL} == 1);
    
    $s_todo .= " $exa $mdl";
    
    #$s_todo .= " > /dev/null" if !$self->{VERBOSE};
    
    #    print "$s_todo\n";
    $self->{CLASSIFY_OUTPUT} = `$s_todo`;
    
    print  $self->{CLASSIFY_OUTPUT} if $self->{VERBOSE};
    
    
    ($self->{CLASSIFY_ERROR}) = $self->{CLASSIFY_OUTPUT} =~ /Leave\-one\-out\ estimate\ of\ the error\:\ error=([\d\.]+)\%/;
    ($self->{CLASSIFY_RECALL}) = $self->{CLASSIFY_OUTPUT} =~ /Leave\-one\-out\ estimate\ of\ the recall\:\ recall=([\d\.]+)\%/;
    ($self->{CLASSIFY_PRECISION}) = $self->{CLASSIFY_OUTPUT} =~ /Leave\-one\-out\ estimate\ of\ the precision\:\ precision=([\d\.]+)\%/;
    
    
    #print "e=$1\n";
    print "$s_todo\n" if $self->{VERBOSE};
    
    
}

sub getError {
    
        my ($self)  = @_;

	return $self->{CLASSIFY_ERROR};

}


sub getRecall {
    
        my ($self)  = @_;

	return $self->{CLASSIFY_RECALL};

}


sub getPrecision {
    
        my ($self)  = @_;

	return $self->{CLASSIFY_PRECISION};

}

sub setKernel {
    my ($self, $k)  = @_;
        
    $self->{K} = $k;
    
}




sub setGamma {
    my ($self, $g)  = @_;
        
    $self->{G} = $g;
    
}


sub setC {
    my ($self, $g)  = @_;
        
    $self->{C} = $g;
    
}



sub addExample {
    
    my ($self, $a_ref, $c, $w, $l)  = @_;
	

   
    
    push @{$self->{RAWEXAMPLES}}, $a_ref;

    if (defined($c)) {
	push @{$self->{CLASSES}},  $c;
    }

    push @{$self->{LABELS}}, $l if (defined($l));
    push @{$self->{WEIGHTS}}, $w if (defined($w));
    

}






sub getExample {
    
    my ($self, $i)  = @_;
	
    
    

    return $self->{EXAMPLES}->[$i];
  
}


sub emptyExamples {

    my ($self)  = @_;
	
    $self->{EXAMPLES} = undef;
    $self->{CLASSES}  = undef;
    
}


# r=1 if if all examples are labeled "0" 
sub _writeExamples {

    my ($self)  = @_;

    for (my $j=1; $j<=$self->{CROSSSVALSIZE}; $j++) {
	
	# write set of examples $i to disk 
	
	open OUT, ">/tmp/svmlight.exa.$j";
	
	my $n = scalar(@{$self->{RAWEXAMPLES}});
	my $cnt = 0;
	for (my $i=0; $i<$n; $i++) {

	    

	    if ($self->{AFFECTATIONS}->[$i] == $j) {

		my $e = $self->{RAWEXAMPLES}->[$i];
		my $l = $self->{LABELS}->[$i];
		
		my $c = (defined($self->{CLASSES}->[$i])?$self->{CLASSES}->[$i]:0);

		for (my $k=0; $k<$self->{WEIGHTS}->[$i]; $k++) {

		    
		    
		    print OUT "$c " . $self->_exampleToSVMLightFormat($e);
		    print OUT "\n";

		    # create a array (class, label)
		    my @a_tmp = ($c, $l);
		    push @{ $self->{WRITTEN}->[$j]}, \@a_tmp;

		}
		
		$cnt++;
	    }
	}

	close OUT;

	if ($self->{VERBOSE} == 1) {
	    print "Written $cnt indiv. to file $j ..\n";
	}
    }
    
    
}

#
# 
#
sub WriteToR {

    my ($self, $efile, $pos, $neg)  = @_;

    $pos = (defined($pos)?$pos:"1");
    $neg = (defined($neg)?$neg:"0");
    
    #
    #  output examples
    #
    open OUT, ">$efile";
    
    my $n = scalar(@{$self->{RAWEXAMPLES}});
    print "n=$n\n";
    for (my $i=0; $i<$n; $i++) {
	my $e = $self->{RAWEXAMPLES}->[$i];
	my $c = $self->{CLASSES}    ->[$i];
	#my $l = $self->{LABELS}     ->[$i];

    	
	if ($c == 1) {
	    $c = $pos;
	} else {
	    $c = $neg;
	}
	
	print OUT "$c\t";
	print OUT join("\t", @$e); print OUT "\n";
	
    }
    close OUT;


}


sub getNbVars {

    my ($self) = @_;


    my $n2 = 0;
    
    if (scalar(@{$self->{MASK}}) != 0) {
	
	for(my $i=0; $i<scalar(@{$self->{MASK}}); $i++) {

	    $n2 ++ if ($self->{MASK}->[$i] == 1);
	}

	$n2 =$n2 * 20;
	

    } else {
	$n2 = length($self->{EXAMPLES}->[0]) * 20;
	
    }
    
    

    return $n2;
}

# $a_ref_index contains the index of the examples that should be written
sub WriteToSNNS {

    my ($self, $file, $a_ref_index)  = @_;

    

    open OUT, ">$file";
    
    print OUT "SNNS pattern definition file V3.2\n";
    print OUT "generated at Mon Apr 25 15:58:23 1994\n";

    print OUT "\n\n";

    my $n1 = 0;

    if (scalar(@$a_ref_index) == 0) { 
       $n1 = scalar(@{$self->{EXAMPLES}});
    } else {
       $n1 = scalar(@$a_ref_index);
    }

    # calc the number of input units

    my $n2 = 0;
    
    if (scalar(@{$self->{MASK}}) != 0) {
	
	for(my $i=0; $i<scalar(@{$self->{MASK}}); $i++) {

	    $n2 ++ if ($self->{MASK}->[$i] == 1);
	}

	$n2 =$n2 * 20;
	
    } else {
	$n2 = length($self->{EXAMPLES}->[0]) * 20;
    }

    print OUT "No. of patterns : $n1\n";
    print OUT "No. of input units : $n2\n";
    print OUT "No. of output units : 1\n";

    print OUT "\n\n";

    my $i = 0;
    foreach my $e ( @{$self->{EXAMPLES}} ) {
	
	if ($a_ref_index && !$self->_in_array($i, @$a_ref_index)) {
	    $i++;
	    next;
	    
	}

	my $t = $self->_sparse_coding($e);
	
	my $c = ($self->{CLASSES}->[$i]>0?1:0);
	
	print OUT "# Input pattern " . ($i+1) . ":\n";
	print OUT $t;
	print OUT "\n";

	print OUT "# Output pattern " . ($i+1) . ":\n";
	print OUT $c;
	print OUT "\n";

	$i++;
	
    }
    
    close OUT;

}







sub _analyzeResults {
    my ($self, $f)  = @_;

    
}

sub classify {

    my ($self)  = @_;

    # write all examples to disk
    $self->_writeExamples();

    $self->{PREDICTIONS} = [];
	

    # run SVM light

    my $s_todo = "svm_classify /tmp/svmlight.exa " .
	$self->{MODEL_FILE} . " /tmp/svmlight.out";

    
     $s_todo .= " > /dev/null"  if  !$self->{VERBOSE};

    print "$s_todo\n" if  $self->{VERBOSE};
    system(($s_todo)) == 0 or die "Cannnot exec svm_classify ..\n";

    
    # get the results

    open IN, "/tmp/svmlight.out" or die "cannot open result file\n";
    my @a = <IN>;

    chomp @a;
    close IN;
    
    my $i_cnt_neg = 0;
    my $i_cnt_pos = 0;

    my $i = 0;
    foreach my $v (@a) {
	
	
	
	$i_cnt_pos++ if ($v > 0.0);
	$i_cnt_neg++ if ($v < 0.0);

	my @a = ($self->{LABELS}->[$i], $self->{CLASSES}->[$i], $v);

	push @{$self->{PREDICTIONS}}, \@a; 
    
	$i ++;
    }
    
    

    #print "Found $i_cnt_pos positifs, $i_cnt_neg negatives ..\n";

    $self->{NBNEG} = $i_cnt_neg;
    $self->{NBPOS} = $i_cnt_pos;


}

#
#  classify a given example file using a givem model file
#
sub _classify {

    my ($self, $exa, $mdl)  = @_;
  
    my $home = `echo \$HOME`; chomp $home;
    
    my $s_todo = "$home/PERL_MODULES/PROGRAMS/SVMLIGHT/svm_classify $exa $mdl /tmp/svmlight.out";

    $s_todo .= " > /dev/null"  if  !$self->{VERBOSE};
    
    print "$s_todo\n" if  $self->{VERBOSE};
    system(($s_todo)) == 0 or die "Cannnot exec svm_classify ..\n";

    
    return Sets::readSet("/tmp/svmlight.out");
    
}


sub setVerbose {

     my ($self, $i)  = @_;
    
    
     $self->{VERBOSE}  = $i;
}

sub setMask {
    
    my ($self, $s)  = @_;
    
    my @a = split //, $s;
    
    
    foreach my $m (@a) {
	push @{$self->{MASK}}, $m;
    }
    
    
}


sub setCrossval {

    my ($self, $i)  = @_;
    
    
    $self->{CROSSVAL}  = $i;

    
    
}

sub setModel {

    my ($self, $s)  = @_;
    
    
    $self->{MODEL_FILE}  = $s;

    
    
}

sub _in_array() {
    my $self = shift(@_);
    my $val = shift(@_);

    #print "is $val in @_ ?\n";

    foreach my $elem(@_) {
        if($val == $elem) {
            return 1;
        }
    }
    return 0;
}


sub _exampleToSVMLightFormat {

    my ($self, $a_ref) = @_;
    my $v = "";
    for (my $i=0; $i<scalar(@$a_ref); $i++) {
	$v .= ($i+1) . ":" . $a_ref->[$i] . " ";
    }
    chop $v;
    return $v;
    
    
}

sub _tosvm {

    my ($self, $s) = @_;

    my @a = split /\ /, $s;

    my $v = "";

    for (my $i=0; $i<scalar(@a); $i++) {
	
	$v .= ($i+1) . ":" . $a[$i] . " ";
	
    } 

    chop $v;

    return $v;
    
}



sub _sparse_coding {

    my ($self, $s) = @_;

my %index = (
    'A' => 0,
    'C' => 1,
    'D' => 2,
    'E' => 3,
    'F' => 4,
    'G' => 5,
    'H' => 6,
    'I' => 7,
    'K' => 8,
    'L' => 9,
    'M' => 10,
    'N' => 11,
    'P' => 12,
    'Q' => 13,
    'R' => 14,
    'S' => 15,
    'T' => 16,
    'V' => 17,
    'W' => 18,
    'Y' => 19);

    my @a = split //, $s;

    my $c = "";


    my $p = 0;
    foreach my $l (@a) {
	
	if ((scalar(@{$self->{MASK}}) > 0) && ($self->{MASK}->[$p] == 0)) {

	    $p++;
	    next;
	}
	#print "l=$l\n";
	
	for(my $i=0; $i<20; $i++) {
	    
	    if ($l eq " ") {
		$c .= "0.0";
	    } elsif ($index{$l} == $i) {
		$c .= "1.0";
	    } else {
		$c .= "0.0";
	    }

	    
	    $c .= " ";
	    
	}
	#print "\n";
	$p++;
    }
    
    chop $c;
    
    return $c;
    
}



1;

