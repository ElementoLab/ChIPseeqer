package ClustalW;
use Sets;
# use strict;


sub new {
    my ($self)     = {};
    $self->{FILE}  = undef;
    $self->{TEXT}  = undef;
    $self->{SEQS}  = undef; #\%tmp;
    $self->{STARS} = "";
    $self->{NAMES} = undef;
    $self->{SEQ_ARRAY} = [];
    
    $self->{NUMCHARNAME}   = undef; 26;
    $self->{SEQ_UA}    = [];
    $self->{NAM_UA}    = [];

    bless($self);
    return $self;
}

sub setNumCharName {
  my ($self, $n) = @_;
  $self->{NUMCHARNAME}   = $n;  
}


sub setSequences {
  my ($self, $a_ref_n, $a_ref_s) = @_;

  $self->{NAMES}     = $a_ref_n;
  $self->{SEQ_ARRAY} = $a_ref_s;
  for (my $i=0; $i<@{$self->{NAMES}}; $i++) {
    $self->{SEQS}->{ $a_ref_n->[$i] } = $a_ref_s->[$i];
  }

}


sub setSequencesToAlign {
  my ($self, $a_ref_n, $a_ref_s) = @_;

  $self->{NAM_UA}  = $a_ref_n;
  $self->{SEQ_UA} = $a_ref_s;
}


sub removeTrailingSequences {

  my ($self, $n, $where) = @_;

  for (my $i=0; $i<@{$self->{NAMES}}; $i++) {

    if ($where == 0) {
      
      #print "BEF ($n):\n";
      #print $self->{SEQS}->{$self->{NAMES}->[$i]};
      #print "\nAFT\n";

      $self->{SEQS}->{$self->{NAMES}->[$i]} =~ s/^.{$n}//;
      
      #print $self->{SEQS}->{$self->{NAMES}->[$i]};
      #print "\n";
      
    } else {
      $self->{SEQS}->{$self->{NAMES}->[$i]} =~ s/.{$n}$//;
    }

  }

}


sub delSeq {
  my ($self, $n) = @_;

  my @a_t = ();
  
  foreach my $r (@{$self->{NAMES}}) {
    if ($r ne $n) {
      push @a_t, $r;
    } 
  }

  $self->{NAMES} = \@a_t;  
  $self->{SEQS}->{$n} = undef;

}


sub getClustalWformat {
  
  my ($self) = @_;
  my $ww  = 80;

  #my $l   = length($self->{SEQ_ARRAY}->[0]);
  my $l   = length($self->{SEQS}->{ $self->{NAMES}->[0] } );

  my $n   = scalar(@{ $self->{NAMES} });
  my $txt = "CLUSTAL W(1.4) multiple sequence alignment\n";
  $txt   .= "\n\n";

  my $nb_blocks = int($l / $ww) + 1;

  for (my $i=0; $i<$nb_blocks; $i++) {
    
    for (my $j=0; $j<$n; $j++) {
      
      $txt .= $self->{NAMES}->[$j];
      $txt .=  substr("                                   ", 0, $self->{NUMCHARNAME} - length($self->{NAMES}->[$j]));
      
      #$txt .= substr($self->{SEQ_ARRAY}->[$j], $i*$ww, $ww); 
      $txt .= substr($self->{SEQS}->{ $self->{NAMES}->[$j] }, $i*$ww, $ww); 
      
      #for (my $k=0; $k<5; $k++) {
      #	$txt .= substr($self->{SEQ_ARRAY}->[$j], $i*$ww + $k*10, 10); 
      #	#  $txt .= " " if ($k != 4);
      #}

      $txt .= "\n";
      

    }
    $txt .= "\n\n";
  }
  
  return $txt;
  
}


sub getPhylipFormat {
  
  my ($self, $h_ref_code) = @_;


  my $l   = length($self->{SEQS}->{ $self->{NAMES}->[0] } );

  my $n   = scalar(@{ $self->{NAMES} });

  printf("%6d%6d\n", $n, $l);


  my $cnt = 10000;
  my $txt = "";
  for (my $j=0; $j<$n; $j++) {

    my $name = $self->{NAMES}->[$j];
    if (!defined($h_ref_code)) {
      $name = substr($name, 0, 10);
    } else {
      my $tmpname = $name;

      $name = $cnt;
      $name .= "-T1" if ($tmpname =~ /T1$/);
      $name .= "-T2" if ($tmpname =~ /T2$/);
      
      $h_ref_code->{$name} = $tmpname;
      $cnt++;
    }

    $txt .= $name;
    $txt .=  substr("                                   ", 0, 10 - length($name));
    $txt .= $self->{SEQS}->{ $self->{NAMES}->[$j] };
    $txt .= "\n";
  
  }
  
  
  return $txt;
  
}



sub getFastaFormat {
  
  my ($self) = @_;
  
  my $txt = "";
  foreach my $n (@{$self->{NAMES}}) {
    my $ss = $self->{SEQS}->{$n};
    $ss =~ s/\-//g;
    $txt .= ">$n\n$ss\n\n";
  }
  return $txt;
}


sub setFile {
    my ($self, $s) = @_;
    $self->{FILE} = $s;
    $self->{HANDLE} = Sets::getRandomString("HANDLE");

    open  $self->{HANDLE}, $self->{FILE} or die "Cannot open $self->{FILE} ..\n";

    #system("cat $self->{FILE}\n");

    my %tmp = ();
    $self->{SEQS} = \%tmp;
   
    $self->{NAMES} = [];
    $self->{STARS} = "";
    $self->{TEXT} = undef;

    $self->_process;   

    close $self->{HANDLE};
}

sub setText {
    my ($self, $s) = @_;
    $self->{TEXT} = $s;
    
    

    my @a = split /[\r\n]/, $s;
    $self->{TEXT_ARRAY} = \@a;
    $self->{CNT} = 0;
    
    $self->_process;   
    
    
}

sub _readline {
    my ($self) = @_;
    
    if ($self->{FILE}) {
	my $tmp = $self->{HANDLE};
	return <$tmp>;
    } else {
	return $self->{TEXT_ARRAY}->[ $self->{CNT} ++];
    }
    
}


sub run {
  my ($self) = @_;
  
  if (defined($self->{SEQ_UA})) {
    
    # create a temp file
    my $tmpfile = Sets::getTempFile("/tmp/toto");
    $tmpfile .= ".seq";

    open OUTCL, ">$tmpfile";
    my $i = 0;
    foreach my $s (@{ $self->{SEQ_UA} }) {
      my $n = $self->{NAM_UA}->[$i];
      print OUTCL ">$n\n$s\n\n";
      $i++;
    }
    close OUTCL;

    my $todo = "clustalw $tmpfile -OUTORDER=INPUT > /dev/null";
    system($todo);
    
    my $outfile = $tmpfile;
    $outfile =~ s/\.seq$/\.aln/;
    
    $self->setFile($outfile);
    
  } else {
    die "no sequences to align\n";
  }
  
}

sub _process {
  my ($self) = @_;

  $self->{STARS} = "";
  my $cnt = 0;
  while (defined(my $l = $self->_readline)) {
    chomp $l;
    next if ($l =~ /CLUSTAL/);
    next if ($l eq "");

    # 16
    my $nn = undef;
    if (defined($self->{NUMCHARNAME})) {
      $nn = $self->{NUMCHARNAME};
    } else {
      my ($side) = $l =~ /^([^\ ]+\ *)/;
      #print "side='$side'\n";
      $nn = length($side);
      $self->{NUMCHARNAME} = $nn;
    }

    my ($n, $s) = $l =~ /^(.{$nn})(.+)$/;

    $n =~ s/\s//g;
    $s =~ s/\s//g;
    
    if ($n eq "") {
      $self->{STARS} .= $s;
    } else {
      if (!defined($self->{SEQS}->{$n})) {
	$self->{NAMES}->[ $cnt ] = $n;
	$cnt ++;
      }
      
      $self->{SEQS}->{$n} .= $s;
      
    }
  
  }

  
  #
  # transform hash into array
  #
  
  my @a_tmp = ();
  for (my $i=0; $i<@{$self->{NAMES}}; $i++) {
    push @a_tmp, $self->{SEQS}->{$self->{NAMES}->[$i]};
  }

  $self->{SEQ_ARRAY} = \@a_tmp;

}

sub lc {
    my ($self) = @_;

    

    foreach my $k (keys(%{ $self->{SEQS} })) {
      $self->{SEQS}->{$k} = lc ($self->{SEQS}->{$k});
    }
	
    
}

sub getHTMLAlignment {
    my ($self) = @_;

    my $i_linewidth = 90;
    my $i = 0;
    my $s_output = "";

    my @a = keys(%{ $self->{SEQS} });


    my $l = length($self->{SEQS}->{$a[0]});
    

    $s_output .= "<table>\n";
    while ($i < $l) {

        my $i_toend = $l - $i;


        foreach my $k (keys(%{ $self->{SEQS} })) {
            $s_output .= "<tr>\n";
            $s_output .= "<td width=\"50\"><tt>\n";
            $s_output .= $k;
            $s_output .= "</tt></td>\n";
            $s_output .= "<td><tt><table  border=0 cellspacing=0 cellpadding=0><tr><td>";
	    
	    my $subseq = substr($self->{SEQS}->{$k}, $i, ($i_linewidth<$i_toend?$i_linewidth:$i_toend));
	    
	    if ($self->{REDIFY} == 1) {
		$subseq =~ s/([ATGC]{1,})/\<\/td\>\<td bgcolor\=\"red\"\>\<font color=\"white\"\>$1\<\/font\>\<\/td\>\<td\>/g;	
	    }

            $s_output .= $subseq;
            $s_output .= "</td></tr></table></tt></td>\n";
            $s_output .= "\n";
            $s_output .= "</tr>\n";
        }

        $s_output .= "<tr><td width=\"50\">&nbsp;</td><td><tt>";
        my $s_colorized = substr($self->{STARS}, $i, ($i_linewidth<$i_toend?$i_linewidth:$i_toend));
	$s_colorized =~ s/\ /\&nbsp\;/g;
        $s_output .= $s_colorized;
        $s_output .= "</tt></td></tr>\n";

        $s_output .= "\n";

        $s_output .= "\n";

        $i += $i_linewidth;

    }
   
     $s_output .= "</table>\n";
    
    
    return $s_output;
    
}

sub kmerize {
    my ($self, $a_ref_kmers) = @_;

    foreach my $k (keys(%{ $self->{SEQS} })) {
	
	foreach my $r (@$a_ref_kmers) {
	    my $se              =  $r->[0];
	    my $suc             =  uc($se);
	    $self->{SEQS}->{$k} =~ s/$suc/$se/ig;
	    my $se              =  Sets::getComplement($se);
	    my $suc             =  uc($se);
	    $self->{SEQS}->{$k} =~ s/$suc/$se/ig;
	} 

    }
}

sub redify {
    my ($self) = @_;

    $self->{REDIFY} = 1;

}
   
    


sub getSeqArray {
    my ($self) = @_;
    
    my @a = ();
    foreach my $k (keys(%{$self->{SEQS}})) {
	push @a, $self->{SEQS}->{$k};
    }

    return \@a; #$self->{SEQS};
}


sub getInputOrderedSeqArray {
    my ($self) = @_;
        
    my @a = ();
    foreach my $k (@{ $self->{NAMES} }) {
	push @a, $self->{SEQS}->{$k};
    }
    
    return \@a; 
}

sub getSeqsWithNames {
  my ($self) = @_;
  
  my @a = ();
  foreach my $k (@{ $self->{NAMES} }) {
    my @tmp = ($k, $self->{SEQS}->{$k});	
    push @a, \@tmp;
  }
  
  return \@a; 
  
}

sub getSeqs {
    my ($self) = @_;
    
    return $self->{SEQS};
    
}


sub getStars {
    my ($self) = @_;
    
    return $self->{STARS};
}



1;
