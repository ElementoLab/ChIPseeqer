package Table;
#use Sets;



sub new {
    my ($class) = @_;

    my ($self) = {};

    $self->{FILE}   = undef;
    $self->{ARRAY}  = [];
    $self->{HANDLE} = undef;

    $self->{LIMIT}  = undef;
    $self->{DELIM}  = "\t";    
    $self->{UC}     = undef;
    $self->{SKIPCOMMENTS} = 0;
    $self->{PHENO}  = {};
    $self->{PHENO_SAMPLE_ORDER} = [];
    $self->{HEADER}   = undef;
    $self->{ROWNAMES} = undef;
    bless $self;

    
    
    return $self;

}


sub processHeader {
  my ($self) = @_;

  $self->{HEADER} = CORE::shift @{ $self->{ARRAY} };
  
  #return $self->{HEADER};
}


sub shiftHeader {
  my ($self) = @_;

  return CORE::shift @{ $self->{ARRAY} };
}




sub getHeader {
  my ($self) = @_;
  return $self->{HEADER};
}

sub getRowNames {
  my ($self) = @_;
  return $self->{ROWNAMES};
}

sub setPhenotypeTable {
  my ($self, $file) = @_;

  open INP, $file or die "Cannot open phenotype table.\n";
  while (my $l = <INP>) {
    chomp $l;
    my @a = split /\t/, $l, -1;
    $self->{PHENO}->{$a[0]} = $a[1];
    push @{ $self->{PHENO_SAMPLE_ORDER} }, $a[0];
  }
  close INP;

}


sub getPhenotypeColIndices {
  my ($self, $p) = @_;

  die "No header set yet.\n" if (!defined($self->{HEADER}));
  
  my @cols = ();
  for (my $i=0; $i<@{$self->{HEADER}}; $i++) {
    #print $self->{HEADER}->[$i] . "\t" . $self->{PHENO}->{ $self->{HEADER}->[$i] } . "; =~ /$p/ ?\n";

    if ($self->{PHENO}->{ $self->{HEADER}->[$i] } =~ /$p/) {
      #print "YES\n";
      push @cols, $i;
    }
  }

  return \@cols;
}

sub getValuesForPhenotype {
  my ($self, $r, $p) = @_;

  my $a_ref_c = $self->getPhenotypeColIndices($p);
  my @v = ();
  foreach my $i (@$a_ref_c) {
    push @v, $r->[$i];
  }
  return \@v;

}


sub reorderColumnsUsingPhenotype {

  my ($self, $a_ref_p_order) = @_;

  my $a_ref_tmp = [];  # temporary matrix

  if ($self->{VERBOSE} == 1) {
    print "Reorganizing matrix cols using coldesc.. \n";
  }

  my $ng     = scalar( @{$self->{ARRAY}} );
  my $nc     = scalar( @{$self->{ARRAY}->[0] } );
 

  # make the order
  #  1. get all labels, assign a unique id to each label
  my %H      = ();
  my $cnt    = 0;
  my @labels = values( %{ $self->{PHENO} } );
  foreach my $l (@labels) {
    if (!defined($H{$l})) {
      $H{ $l } = $cnt;
      #print "$l\n";
      $cnt++;
    }
  }	
  
  my @order = ();
  $order[0] = 0;
  
  my @labels_order = undef;
  if (!defined($a_ref_p_order)) {
    @labels_order = sort {$b cmp $a}  keys(%H);
  } else {
    @labels_order = @$a_ref_p_order;
  }

  foreach my $lab (@labels_order) {
    for (my $i=1; $i<@{$self->{HEADER}}; $i++) {
      if ( $self->{PHENO}->{ $self->{HEADER}->[$i] } eq $lab) {
	push @order, $i;
      }
    }
  }

  #print "Order has " . scalar(@order) . " elements\n";
  
      

  for (my $i=0; $i<$ng; $i++) {

    foreach my $j (@order) {
      push @{$a_ref_tmp->[$i]}, $self->{ARRAY}->[$i]->[$j];
    }

    #for (my $j=0; $j<$nc; $j++) {
    #  my $newj = $order[$j];      
    #  $a_ref_tmp->[$i]->[$j] = $self->{ARRAY}->[$i]->[$newj];
    #}
    
  }  
  
  $self->{ARRAY} = $a_ref_tmp;
  
  my $tmp = [];

  foreach my $j (@order) {
    push @$tmp,  $self->{HEADER}->[$j] . " " . $self->{PHENO}->{ $self->{HEADER}->[$j] };
  }

  #for (my $j=0; $j<$nc; $j++) {
  #  my $newj   = $order[$j];
  #  $tmp->[$j] = $self->{HEADER}->[$newj] . " " . $self->{PHENO}->{ $self->{HEADER}->[$newj] };
  #}
  $self->{HEADER} = $tmp;

}


sub setSkipComments {
  my ($self, $d) = @_;
  $self->{SKIPCOMMENTS} = $d;
}



sub setDelim {
	my ($self, $d) = @_;
	$self->{DELIM} = $d;
}

sub shift {
    my ($self) = @_;
    return CORE::shift(@{ $self->{ARRAY} });
}

#
# get the number of rwos
#
sub getNbRows {
    my ($self) = @_;

    return scalar( @{ $self->{ARRAY} } );

}


sub setUC {
    my ($self, $l) = @_;
    $self->{UC} = $l;
}


sub setLimit {
    my ($self, $l) = @_;
    $self->{LIMIT} = $l;
}



sub setFile {
    my ($self, $file) = @_;

    

    $self->{FILE} = $file;
    
    #print "self->{FILE}=$self->{FILE}\n";

}


sub nextRow {
     my ($self) = @_;

    if (!defined($self->{HANDLE})) {
	$self->{HANDLE} = Sets::getRandomString("HANDLE");
	
	die "nextRow: Please define a file ..\n" if (!defined($self->{FILE}));
					    
	open $self->{HANDLE}, $self->{FILE} or die "Cannot open $self->{FILE}\n";


    }

    
    my $l = $self->_readline;

    if ($l) {
	chomp $l;
        my $d = $self->{DELIM};
	my @a = split /$d/, $l;
	
	return \@a;
    } else {

	#$self->dispose;
	return undef;
    }
    
}


sub dispose {
    my ($self) = @_;
    
    if (defined($self->{HANDLE})) {
	close $self->{HANDLE};
    }
    $self->{HANDLE} = undef;
}

sub _readline {
    my ($self) = @_;
    
    if ($self->{FILE}) {
	my $tmp = $self->{HANDLE};
	return <$tmp>;
    } 
    else {
	die "Please define a file ..\n";
    }
    
}


sub getFastArray {
    my ($file) = @_;

    open IN, $file or die "Cannot open file \"$file\" ..\n";
   
    my @out = ();
    while (my $l = <IN>) {
	chomp $l;
	my @a = split(/\t/, $l, -1);
	push @out, \@a;
    }
    close IN;
    
    return \@out;
}

#
#  load a file containing the table
#
sub loadFile {
    
    my ($self, $file, $rn, $cn) = @_;

    $self->{ARRAY} = [];

    if (!defined($file)) {
	# "Please provide a file name for Table.pm ..\n";
	return 0;
    }

    $self->{FILE} = $file;

    open INO4, $self->{FILE} or  die "Table.pm: cannot open file \"$self->{FILE}\" ..\n";

    if (defined($cn) && ($cn == 1)) {
      #print "Getting header.\n";
      my $l = <INO4>;
      chomp $l;
      my @a = split /\t/, $l, -1;
      if (defined($rn) && ($rn == 1)) {
	CORE::shift @a;
      }
      $self->{HEADER} = \@a;
    }

    if (defined($rn) && ($rn == 1)){
      $self->{ROWNAMES} = [];
    }

    my $cnt = 0;
    while (my $l = <INO4>) {
	chomp $l;

	if (($self->{SKIPCOMMENTS} == 1) && ($l =~ /^\#/)) {
	  next;
	}

	my $d = $self->{DELIM};
	my @a = split(/$d/, $l, -1);

	if ($self->{UC} && ($self->{UC} == 1)) {
	    
	    @a = map(uc, @a); 
	}
	
	if (defined($rn) && ($rn == 1)) {
	  my $cname = CORE::shift @a;
	  push @{$self->{ROWNAMES}}, $cname;
	}

	push @{ $self->{ARRAY} }, \@a;
	
	$cnt ++;

	last if (defined($self->{LIMIT}) && ($cnt == $self->{LIMIT}));
	
    }
    close INO4;

    #print "Array  contains " . scalar(@ {$self->{ARRAY}} ) . "\n"; 
   
    return 1;
}

#
#  get HTML table
#
sub getHtmlTable {
    my ($self) = @_;

    my $txt = "<table border=1>\n";
    
    foreach my $r (@{ $self->{ARRAY} }) {
	$txt .= "<tr>";
	map { $_ = "<td><tt>$_</tt></td>" } @$r;
	my $tmp .= join("", @$r);
	
	$tmp =~ s/\ /\&nbsp\;/g;
	
	$txt .= $tmp;
	
	$txt .= "</tr>\n";
	
    }

    $txt   .= "</table>\n";

    
    
    return $txt;
	
}


#
#  get a Latex table
#
sub getLatexTable {

    my ($self) = @_;


    # how many cols ?
    my $n = scalar( @{ $self->{ARRAY}->[0] } );

    my $c = "l" x ($n+1);


    my $txt = '\\' . "begin{tabular}{$c}" . "\n";
    
    foreach my $r (@{ $self->{ARRAY} }) {

	my $tmp .= join(" & ", @$r);
	

	
	$txt .= $tmp;
	
	$txt .= '\\\\' . "\n";
	
    }

    $txt   .= '\end{tabular}' . "\n";
      
    return $txt;
	
}


#
#  
#
sub getLyxTable {

    my ($self) = @_;

    # number of columns
    my $nc = scalar( @{ $self->{ARRAY}->[0] } );
    
    my $nr = scalar(@{ $self->{ARRAY} });
    

    my $txt = "\\begin_inset  Tabular\n";
    $txt   .= "<lyxtabular version=\"3\" rows=\"$nr\" columns=\"$nc\">\n";
    $txt   .= "<features>\n";
    for (my $i=0; $i<$nc; $i++) {
	$txt .= "<column alignment=\"left\" valignment=\"top\" leftline=\"false\" width=\"0\">\n";
    }
    
    foreach my $r (@{ $self->{ARRAY} }) {
	$txt .= "<row topline=\"true\">\n";
	
	
	foreach my $b (@$r) {
	    $txt .= "<cell alignment=\"left\" valignment=\"top\" topline=\"true\" leftline=\"false\" usebox=\"none\">\n";
	    $txt .= "\\begin_inset Text\n\n";
	    $txt .= "\\layout Standard\n\n";

	    $txt .= "$b\n";
	    $txt .= "\\end_inset\n\n"; 
	    $txt .= "</cell>\n\n";
	}

	$txt .= "</row>\n";
	
    }


    $txt .= "</lyxtabular>\n";

      
    return $txt;

}





#
#   get index by the ith column
#
sub getIndex {
    my ($self, $i) = @_;

    my %hash = ();
    
    foreach my $r (@{ $self->{ARRAY} }) {
	$hash{ $r->[$i] } = $r;
    }
    
    return \%hash;

}


#
#   get index by the ith column
#
sub getIndexShifted {
    my ($self) = @_;

    my %hash = ();
    
    foreach my $r (@{ $self->{ARRAY} }) {
      my @a = @$r;
      my $n = CORE::shift @a;
      $hash{ $n } = \@a;
    }
    
    return \%hash;

}


#
#   get index by the ith column
#
sub getHash {
    my ($self, $i) = @_;

    my %hash = ();
    
    foreach my $r (@{ $self->{ARRAY} }) {
	$hash{ $r->[$i] } = $r;
    }
    
    return \%hash;

}


#
#   get index of $jth column by the ith column
#
sub getIndexColumnsKV {
    my ($self, $i, $j) = @_;

    my %hash = ();
    
    foreach my $r (@{ $self->{ARRAY} }) {
	$hash{ $r->[$i] } = $r->[$j];
    }
    
    return \%hash;

}

#
#   get index of $jth column by the ith column
#
sub getIndexKV {
    my ($self, $i, $j) = @_;

    my %hash = ();
    
    foreach my $r (@{ $self->{ARRAY} }) {
	$hash{ $r->[$i] } = $r->[$j];
    }
    
    return \%hash;

}


sub getArray {
    my ($self) = @_;

    return $self->{ARRAY};
    
}

#
#  
#
sub getColumn {    
    my ($self, $i) = @_;
    my @a = ();    
    foreach my $r (@{ $self->{ARRAY} }) {
	push @a, $r->[$i];
    }    
    return \@a;
}


#
#  
#
sub getRow {    
    my ($self, $i) = @_;

    return $self->{ARRAY}->[$i];

}


#
#  
#
sub getNbColumns {
    my ($self) = @_;

    return scalar(@{ $self->{ARRAY}->[0] });

}

#
#   get an array of hashes
#
sub getArrayOfHashes {
    my ($self, @a) = @_;
    
    
    my @a_hash = ();
    
    foreach my $r (@{ $self->{ARRAY} }) {
	my %h = ();
	my $i = 0;
	foreach my $e (@$r) {
	    $h{ $a[$i] } = $e;
	    $i++;
	}
	push @a_hash, \%h;
    }
    return \@a_hash;
}

sub sortByCol {
  my ($self, $i, $order) = @_;
  $self->sortbycol($i, $order);
}

sub sortbycol {
  my ($self, $i, $order) = @_;

  if ($order == 1) {
    @{ $self->{ARRAY} } = sort { $a->[$i] <=> $b->[$i] } @{ $self->{ARRAY} };
  } else {
    @{ $self->{ARRAY} } = sort { $b->[$i] <=> $a->[$i] } @{ $self->{ARRAY} };
  }
  
}

sub randomizeColumn {
    my ($self, $n) = @_;

    my $a_ref = $self->getColumn($n);
    my $a_ref_shu = Sets::shuffle_array($a_ref);
    my $nn = scalar(@{ $self->{ARRAY} });
    for (my $i=0; $i<$nn; $i++) {
	$self->{ARRAY}->[$i]->[$n] = $a_ref_shu->[$i];
    }

}


sub printTable {
    my ($a_ref) = @_;

    foreach my $r (@$a_ref) {
	print join("\t", @$r) . "\n";
    }

}


sub print {
    my ($self) = @_;
    
    foreach my $r (@{ $self->{ARRAY} }) {
	print join("\t", @$r) . "\n";
    }

}


sub getBidimensionalHash {
  my ($self, $f) = @_;

  my %H = ();

  open IN, $f or die "can't open $f ..\n";
  my $l = <IN>;
  chomp $l;
  my @a = split /\t/, $l, -1;

  CORE::shift @a;
  while (my $l = <IN> ) {

     chomp $l;
     my @b = split /\t/, $l, -1;
     
     my $n = CORE::shift @b;
  
     my $t = scalar(@b);
     for (my $i=0; $i<$t; $i++) {
       
       $H{ $n } { $a[$i] } = $b[$i];
       
     }

  }

  close IN;
  
  return \%H;

}

1;
