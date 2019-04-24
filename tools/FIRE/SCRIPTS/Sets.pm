package Sets;
#use strict;

use Table;
#
# get a distribution from a set by binning
#    array of (middle, frequency)
#

srand;

sub new {
    my $self  = {};
    bless($self);
    return $self;
}


sub transpose {
  my ($a_ref_m) = @_;
  
  my $n = @$a_ref_m;
  my $m = @{$a_ref_m->[0]};

  my $a_ref_m2 = [];
  for (my $i=0; $i<$n; $i++) {
    for (my $j=0; $j<$m; $j++) {
      $a_ref_m2->[$j][$i] = $a_ref_m->[$i][$j];
    }
  }
  return $a_ref_m2;

}


sub getpwd {  
  my $pwd = `pwd`; 
  chomp $pwd;
  return $pwd;
}


sub count_file_lines {
  my ($f) = @_;

  die "File does not exist\n" if (! -e $f);

  my $out = `wc -l $f`;
  chomp $out;
  $out =~ s/\ .+$//;
  return $out;
}

#
# Rosen's Algorithm 1 for generating all permutations, page 138 (my copy)
#
sub get_all_permutations {
  my ($n) = @_;

  
  my @OUT = ();

  my @a = ();
  for (my $i=0; $i<$n; $i++) {
    push @a, $i;
  }
  my @e = @a; push @OUT, \@e;
    

  my @b = reverse @a; 
  my $txt_b = join("-", @b);
  my $txt_a = "";

  while ($txt_a ne $txt_b) {
    my $m = undef;

    #print "Current perm = " . join(" ", @a) . "\n";

    #
    # found m the rightmost location such that a[m] is smaller than a[m+1]
    #
    for ($m=$n-2; $m>=0; $m--) {
      #print  $a[$m+1] . " > " . $a[$m] . " ?\n";
      if ($a[$m+1] > $a[$m]) {
	last;  
      } 
    }

    #print " -> m=$m\n";

    # find the smallest value larger than a[m] to the right of a[m]
    my $idx = undef;

    for (my $i=$m+1; $i<$n; $i++) {
      if ($a[$i] > $a[$m]) {
	if (!defined($idx)) {
	  $idx = $i;
	} else {
	  
	  if ($a[$i] < $a[$idx]) {
	    $idx = $i;
	  }

	}
      }
    }

    #print " -> idx=$idx\n";

    my @c = ();
    for (my $i=0; $i<$m; $i++) {
      push @c, $a[$i];
    }
    $c[$m] = $a[$idx];
    
    my @tmp = ();
    for (my $i=$m; $i<$n; $i++) {
      if ($i == $idx) {
	next;
      } else {
	push @tmp, $a[$i];
      }
    }
    @tmp = sort @tmp;
    push @c, @tmp;

    @a = @c;
    
    my @d = @a; push @OUT, \@d;
    
    $txt_a = join("-", @a);
    
  }
  
  return \@OUT;
}



sub swap {
  my ($r1, $r2) = @_;
  my $tmp = $$r1;
  $$r1 = $$r2;
  $rr2 = $tmp;
}

sub nbLinesInFile {
  my ($f) = @_;
  
  open IN, $f;
  my @a = <IN>;
  close IN;
  
  return scalar(@a);
}

sub defined_exists_or_die {
  
  my ($v, $m) = @_;
  
  if (!defined($v)) {
    die "$m";
  }
  if (! -e $v) {
    die "$v does not exists.";
  }
  
}

sub consensus2re {
  
  my ($c) = @_;
  
  my @a    = split //, $c;
  
  my @ch   = ();
  my @d_ch = ();
  
  $d_ch[0]  = "N";
  $d_ch[1]  = "A";
  $d_ch[2]  = "C";
  $d_ch[3]  = "G";
  $d_ch[4]  = "T";
  $d_ch[5]  = "M";
  $d_ch[6]  = "R";
  $d_ch[7]  = "W";
  $d_ch[8]  = "S";
  $d_ch[9]  = "Y";
  $d_ch[10] = "K";
  $d_ch[11] = "V";
  $d_ch[12] = "H";
  $d_ch[13] = "B";
  $d_ch[14] = "D";

  $ch[0]    = ".";
  $ch[1]    = "A";
  $ch[2]    = "C";
  $ch[3]    = "G";
  $ch[4]    = "T";
  $ch[5]    = "[AC]";
  $ch[6]    = "[AG]";
  $ch[7]    = "[AT]";
  $ch[8]    = "[CG]";
  $ch[9]    = "[CT]";
  $ch[10]   = "[GT]";
  $ch[11]   = "[ACG]";
  $ch[12]   = "[ACT]";
  $ch[13]   = "[CGT]";
  $ch[14]   = "[AGT]";

  for (my $i=0; $i<@d_ch; $i++) {
    $H{ $d_ch[$i] } = $ch[$i];
  }
  
  while (@a < 9) {
    @a = ('N', @a, 'N');
  }
    
  my $txt = "";
  foreach my $r (@a) {
    $txt .= $H{$r};
  }
  
  return $txt;
}


sub myre2scanace {
  my ($re) = @_;
  
  my $oo = get_array_from_re($re);
  my $m = @$oo;
  
  my $txt = myre2wm($re);
  my $sta = '*' x $m;
  $txt   .= "$sta\n"; 
  
  return $txt;
}

#
#  aa RE => WM
#
sub myre2wm_aa {

  my ($re)  = @_;

  my @aa    = ('A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V');

  my @a     = split //, $re;  
  my @lines = ();
  my $i     = 0; 
  my $n     = length($re);
  my $oo    = get_array_from_re($re);
  my $m     = @$oo;
  
  while ($i < $n) {
    
    if ($a[$i] eq '.') {

      my $l = "";
      
      foreach my $theaa (@aa) {
	$l .= $theaa x 5;
      }

      push @lines, $l;
    }
    
    elsif ($a[$i] eq '[') {
            
      my @b = ();
      $i++;
      while ($a[$i] ne ']') {
        push @b, $a[$i];
        $i ++; 
      }
      
      
      my $m    = scalar( @b );
      my $cbaa = int(100 / $m);

      my $l    = "";
      foreach my $s (@b) {
	$l .= $s x $cbaa;
      }
      
      if (length($l) < 100) {

	while (length($l) < 100) {
	  $l .= 'N';
	}
	
      }

      push @lines, $l;

    } else {
      my $l = $a[$i] x 100;
      push @lines, $l;
    }
    
    $i ++;
  }	
  
  
  my @M = ();
  my $col = 0;
  foreach my $l (@lines) {
    my @a = split //, $l; 
    for (my $i=0; $i<100; $i++) {
      $M[ $i ] [ $col ] = $a[$i];
    }
    $col ++;
  } 
  
  my $mo = "";
  for ( my $i=0; $i<100; $i++) {
    for (my $j=0; $j<$m; $j++) {
      $mo .= $M[$i][$j]; 
    }
    $mo .= "\n";
  }
  
  return $mo;
}



sub myre2wm {
  my ($re) = @_;

  my @a = split //, $re;
  
  my @lines = ();
  
  my $i = 0; my $n = length($re);

  my $oo = get_array_from_re($re);
  my $m = @$oo;
  
  while ($i < $n) {
    
    
    if ($a[$i] eq '.') {
      my $l = "";
      $l .= 'A' x 25;
      $l .= 'C' x 25;
      $l .= 'G' x 25;
      $l .= 'T' x 25;
      push @lines, $l;
    }
    
    elsif ($a[$i] eq '[') {
      
      
      my @b = ();
      $i++;
      while ($a[$i] ne ']') {
        push @b, $a[$i];
        $i ++; 
      }
      
      
      my $m = scalar( @b );
      if ($m == 2) {
        my $l = "";
        foreach my $s (@b) {
          $l .= $s x 50;
        }
        push @lines, $l;
	
      } elsif ($m == 3) {
	my $l = "";
	$l .= $b[0] x 34;
	$l .= $b[1] x 33;
	$l .= $b[2] x 33;
	push @lines, $l;
      }
      
    } else {
      my $l = $a[$i] x 100;
      push @lines, $l;
    }
    
    $i ++;
  }	
  
  
  my @M = ();
  my $col = 0;
  foreach my $l (@lines) {
    my @a = split //, $l; 
    for (my $i=0; $i<100; $i++) {
      $M[ $i ] [ $col ] = $a[$i];
    }
    $col ++;
  } 
  
  my $mo = "";
  for ( my $i=0; $i<100; $i++) {
    for (my $j=0; $j<$m; $j++) {
      $mo .= $M[$i][$j]; 
    }
    $mo .= "\n";
  }
  
  return $mo;
}


sub interp_from_matlab_colormap {
  my ($r, $a_ref_map, $min, $max) = @_;

  if (@$a_ref_map == 0) {
    die "Color map is empty ... exiting.\n";
  }

  if ($r > $max) {
    $r = $max;
  }
  
  if ($r < $min) {
    $r = $min;
  }
  
  my $n = scalar(@$a_ref_map);

  # $n -> $max
  #  ? -> $r
  #  0 -> $min
  my $i = int(0.5 + $n * ($r - $min ) / ($max - $min));

  if ($i < 0) {
    $i = 0;
  }
  
  if ($i > $n-1) {
    $i = $n-1;
  }

  return ( int(0.5+$a_ref_map->[$i]->[0] * 255), int(0.5+$a_ref_map->[$i]->[1] * 255), int(0.5+$a_ref_map->[$i]->[2] * 255) );
  
}  


sub interp_general {

  my ($r, $cc1, $cc2, $min, $max) = @_;

  if ($r > $max) {
    $r = $max;
  }
  
  if ($r < $min) {
    $r = $min;
  }
  
  
  my ($a1, $b1, $c1) = @$cc1;
  my ($a2, $b2, $c2) = @$cc2;
  
  my $a_ref_A = Sets::getInterpLine($a1, $min, $a2, $max);
  my $a_ref_B = Sets::getInterpLine($b1, $min, $b2, $max);
  my $a_ref_C = Sets::getInterpLine($c1, $min, $c2, $max);
  
  my $ai      = int(0.5 + $a_ref_A->[0] * $r + $a_ref_A->[1]);
  my $bi      = int(0.5 + $a_ref_B->[0] * $r + $a_ref_B->[1]);
  my $ci      = int(0.5 + $a_ref_C->[0] * $r + $a_ref_C->[1]);
  
  return ($ai, $bi, $ci); 
  
  #my $A = ($a1 - $a2) / ($min - $max);
  #my $B = ( - $a1 * $max + $a2 * $min ] / ($min - $max);
  
  

}




sub getInterpLine {

  my ($y1, $x1, $y2, $x2) = @_;

  my $A = ($y1 - $y2) / ($x1 - $x2);
  my $B = ( - $y1 * $x2 + $y2 * $x1 ) / ($x1 - $x2);
  
  return [$A, $B];
  
}

sub order {

  my ($a_ref, $r) = @_;
  my @a = ();
  my $cnt = 0;

  foreach my $s (@$a_ref) {
    my @a_tmp = ($s, $cnt);
    push @a, \@a_tmp;
    $cnt ++;
  }
  if (!defined($r) || ($r == 1)) { 
    @a = sort { $a->[0] <=> $b->[0] } @a;
  } else {
    @a = sort { $b->[0] <=> $a->[0] } @a;
  }
  my @b = ();

  foreach my $s (@a) {
    push @b, $s->[1];
  }
  return \@b;

}


sub ranks {

  my ($a_ref, $r) = @_;
  my @a = ();
  my $idx = 0;
  foreach my $s (@$a_ref) {
    my @a_tmp = ($s, $idx);
    push @a, \@a_tmp;
    $idx ++;   
  }
  if (!defined($r) || ($r == 1)) { 
    @a = sort { $a->[0] <=> $b->[0] } @a;
  } else {
    @a = sort { $b->[0] <=> $a->[0] } @a;
  }

  my $rank = 0;
  my @b = ();
  foreach my $s (@a) {
    if (defined($s->[0]) && ($s->[0] ne "")) {
      $b[$s->[1]] = $rank;
      $rank ++;
    } else {
      $b[$s->[1]] = undef;
    }
  }
  return \@b;

}


sub hash_order {
  my ($h_ref) = @_;

  my @a = ();
  foreach my $k (keys(%$h_ref)) {
    my @a_tmp = ($k, $h_ref->{$k});
    push @a, \@a_tmp;
  }
  
  # sort based on values
  @a = sort { $a->[1] <=> $b->[1] } @a;
  
  

  my @b = ();
  foreach my $r (@a) {
    push @b, $r->[0];
  }
  
  return \@b;
}

sub get_rank_from_array {
  my ($a_ref) = @_;

  my @a = ();
  my $i = 0;
  foreach my $r (@$a_ref) {
    my @a_tmp = ($i, $r);
    push @a, \@a_tmp;
    $i ++;
  }
  
  # sort based on values
  @a = sort { $a->[1] <=> $b->[1] } @a;
  
  
  $i = 0;
  my @b = ();
  foreach my $r (@a) {
    my $idx = $r->[0];
    $b[ $idx ] = $i;
    $i ++;
  }
  
  return \@b;
}


sub binom_test_greater {
  
  my ($n, $k, $p, $p1) = @_;

  my $out = `$home/PROGRAMS/MIMOTIFS/binom_test_greater $n $k $p 1`;
  #print "/home/elemento/PROGRAMS/MIMOTIFS/binom_test_greater $n $k $p 1\n";
  $out =~ s/[\n\r]//g;

  my @a = split /\t/, $out, -1;
  
  $out = $a[0];

  $$p1 = $a[1];
  
  return $out;

}


sub getNiceDateTime {
  
  my ($us) = @_;
  
  my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
  my @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
  my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
  my $year = 1900 + $yearOffset;
  my $theTime = "$hour" . "h:$minute" . "m:$second" . "s, $weekDays[$dayOfWeek] $months[$month] $dayOfMonth, $year";
  
  if (defined($us)) {
    $theTime =~ s/\ /\_/g;
    $theTime =~ s/\,//g;
    $theTime =~ s/\:/\_/g;
  }
  
  return $theTime;

}


sub get_array_from_re {
  my ($re) = @_;

  my @a = split //, $re;
  my $i = 0; my $n = length($re);
  my @d = ();
  while ($i < $n) {
   
    if ($a[$i] eq '[') {
	
      my $b = "[";
      $i++;
      while ($a[$i] ne ']') {
        $b .= $a[$i];
        $i ++; 
      }
      $b .= "]";
      push @d, $b;
    } else {
      push @d, $a[$i];
    }
  	 	
  $i ++;
 }	

  return \@d;
}



sub shuffle_re {
  my ($re) = @_;

  my @a = split //, $re;
  my $i = 0; my $n = length($re);
  my @d = ();
  while ($i < $n) {
   
    if ($a[$i] eq '[') {
	
      my $b = "[";
      $i++;
      while ($a[$i] ne ']') {
        $b .= $a[$i];
        $i ++; 
      }
      $b .= "]";
      push @d, $b;
    } else {
      push @d, $a[$i];
    }
  	 	
  $i ++;
 }	

  return shuffle_array(\@d);
}



sub line_with_index {
    my $data_file   = shift;
    my $index_file  = shift;
    my $line_number = shift;

    my $size;               # size of an index entry
    my $i_offset;           # offset into the index of the entry
    my $entry;              # index entry
    my $d_offset;           # offset into the data file

    $size = length(pack("N", 0));
    $i_offset = $size * $line_number;
    seek($index_file, $i_offset, 0) or return;
    read($index_file, $entry, $size);
    $d_offset = unpack("N", $entry);
    seek($data_file, $d_offset, 0);
    return scalar(<$data_file>);
}

# usage:
#$file = $ARGV[0];
#open(FILE, "< $file")         or die "Can't open $file for reading: $!\n";
#open(INDEX, "< $file.idx")
#        or die "Can't open $file.idx for read/write: $!\n";
#$line = line_with_index(*FILE, *INDEX, $ARGV[1]);
#print $line;




sub get_parameter {

  my ($a_ref_params, $param) = @_;
  
  my $argc = scalar(@$a_ref_params);

  my $i = 0;

  while (($i < $argc) && ($param ne $a_ref_params->[$i])) {
    $i++;
  }
  
  if ($i<$argc) {
    return $a_ref_params->[$i+1];
  } else {
    return undef;
  }
}


sub exist_parameter {

  my ($a_ref_params, $param) = @_;
  
  my $argc = scalar(@$a_ref_params);

  my $i = 0;

  while (($i < $argc) && ($param ne $a_ref_params->[$i])) {
    $i++;
  }
  
  if ($i<$argc) {
    return 1;
  } else {
    return 0;
  }
}




#
#  returns connected components out of an adjacency list H{ node } => node
#
sub get_connected_components {
  my ($h_ref_ADJ) = @_;
  
  my @COMP = ();

  # init colors
  my %COL = ();
  foreach my $n (keys(%$h_ref_ADJ)) {
    $COL{$n} = 0;
  }
  
  # DFS
  foreach my $n (keys(%$h_ref_ADJ)) {
    if ($COL{$n} == 0) {
      my @a_nodes = ();
      DFS($h_ref_ADJ, $n, \%COL, \@a_nodes);
      push @COMP, \@a_nodes;
    }
  }

  return \@COMP;
} 
	

sub DFS {
    my ($h_ref_adj, $node, $h_ref_col, $a_ref_nodes) = @_;

    $h_ref_col->{ $node } = 1;
    
    foreach my $next (@{ $h_ref_adj->{ $node } }) {
	next if (!defined($h_ref_adj->{$next}));
	if ($h_ref_col->{ $next } == 0) {
            DFS($h_ref_adj, $next, $h_ref_col, $a_ref_nodes);
        }
    }
    $h_ref_col->{ $node } = 2;
    push @$a_ref_nodes, $node;
}





sub getNbGenesInFastaFile {
  my ($f) = @_;

  open INK, $f or die "cannot open $f\n";
  my @a = <INK>;
  my @b = grep /^\>/, @a;
  close INK;
  return scalar(@b);
}


sub getDictionnary {
  my ($f) = @_;
  my %H = ();
  open IN, $f;
  while (my $l = <IN>) {
    chomp $l;
    my @a = split /\t/, $l;
    
    my $n = shift @a;

    $H{ $n } = \@a;
    
    
  }
  close IN;
  
  return \%H;
}

sub loadGeneTable {
    
    my ($file, $h_ref_CHR, $h_ref_BOUNDARIES_P, $h_ref_STRAND, $h_ref_BOUNDARIES_T) = @_;
    
    my $ta = Table->new;
    $ta->loadFile($file);
    my $a_ref = $ta->getArrayOfHashes( ("ORF", "SCAFFOLD", "START_P", "END_P", "STRAND", "START_T", "END_T") );
        
    foreach my $r (@$a_ref) {

	if ($r->{"END_T"} < $r->{"START_T"}) {
	    my $tt = $r->{"START_T"};
	    $r->{"START_T"} = $r->{"END_T"};
	    $r->{"END_T"} = $tt;
	}
	
	if ($r->{"END_P"} < $r->{"START_P"}) {
	    my $tt = $r->{"START_P"};
	    $r->{"START_P"} = $r->{"END_P"};
	    $r->{"END_P"} = $tt;
	}

	if (!defined($r->{"START_T"})) {
	    $r->{"START_T"} = $r->{"START_P"};
	    $r->{"END_T"}   = $r->{"END_P"};
	}

	# chromosome
	$h_ref_CHR->{ $r->{"ORF"} }             = $r->{"SCAFFOLD"}; 

	# protein
	$h_ref_BOUNDARIES_P->{ $r->{"ORF"} }->[0] = (defined($h_ref_BOUNDARIES_P->{ $r->{"ORF"} }->[0])?Sets::min($h_ref_BOUNDARIES_P->{ $r->{"ORF"} }->[0], $r->{"START_P"}):$r->{"START_P"});
	$h_ref_BOUNDARIES_P->{ $r->{"ORF"} }->[1] = (defined($h_ref_BOUNDARIES_P->{ $r->{"ORF"} }->[1])?Sets::max($h_ref_BOUNDARIES_P->{ $r->{"ORF"} }->[1], $r->{"END_P"}):$r->{"END_P"});
	
	# strand
	$h_ref_STRAND->{ $r->{"ORF"} }          = $r->{"STRAND"};
	
	# transcript
	$h_ref_BOUNDARIES_T->{ $r->{"ORF"} }->[0] = (defined($h_ref_BOUNDARIES_T->{ $r->{"ORF"} }->[0])?Sets::min($h_ref_BOUNDARIES_T->{ $r->{"ORF"} }->[0], $r->{"START_T"}):$r->{"START_T"});
	$h_ref_BOUNDARIES_T->{ $r->{"ORF"} }->[1] = (defined($h_ref_BOUNDARIES_T->{ $r->{"ORF"} }->[1])?Sets::max($h_ref_BOUNDARIES_T->{ $r->{"ORF"} }->[1], $r->{"END_T"}):$r->{"END_T"});
    }

}




sub translateSet {
    my ($set, $dic) = @_;

    my @a = ();
    foreach my $r (@$set) {
	if (defined($dic->{$r})) {
	    push @a, $dic->{$r};
	}
    }
    return \@a;
}



 

sub isMiRNATarget {
    my ($f, $sp) = @_;
    
    my $file = undef;
    if ($sp eq "fly") {
	$file = "/home/olly/DATA/DROSOPHILA/MIRNA/5p_kmers_miRNAs.txt";
    }

    my %H = ();
    open IN, $file or die "no mirna file\n";
    while (my $l = <IN>) {
	chomp $l;
	my @a = split /\t/, $l, -1;
	my $n = shift @a;
	$H{ $n } = \@a;
    }
    close IN;

    return $H{ $f }; 

}


sub unlink_blast_files {

    my ($f) = @_;

    unlink "$f.nhr" if (-e "$f.nhr");
    unlink "$f.nin" if (-e "$f.nin");
    unlink "$f.nsd" if (-e "$f.nsd");
    unlink "$f.nsi" if (-e "$f.nsi");
    unlink "$f.nsq" if (-e "$f.nsq");


}

sub gccontent {
    my ($s) = @_;

    my @a = split //, $s;

    my %H = ();
    foreach my $n (@a) {
	$H{ $n } ++;
    }

    my $gc = ($H{G} + $H{C}) / ($H{A} + $H{T} + $H{G} + $H{C});

    return $gc;
    
    
}


sub initHash {
    my %h = ();

    return \%h;
}

sub cmdLine {
    my ($argcmin, $s) = @_;
    our @ARGV;

    #print scalar(@ARGV);
    
    

    if (scalar(@main::ARGV) < $argcmin) {
	die "Usage : $^X $s\n";
    }
    
    
}

#
#  return the maximum overlap between 2 sequences
#
sub seqOverlap {
    my ($s1, $s2, $ref) = @_;

    my $l1 = length($s1);
    my $l2 = length($s2);
    my $best_ov = 0;
    for (my $i=1; $i<=$l1; $i++) {
	my $ss1 = substr($s1, 0, $i);
	$ss1 .= '$';
	#print " $ss1\n";
	if ($s2 =~ /$ss1/) {
	    if ($i > $best_ov) {
		$best_ov = $i;
		$$ref = -1;
	    }
	}

	$ss1 = substr($s1, $l1-$i, $i);
	$ss1 = '^' . $ss1;
	#print " $ss1\n";
	if ($s2 =~ /$ss1/) {
	    if ($i > $best_ov) {
		$best_ov = $i;
		$$ref = +1;
	    }
	    #$best_ov = $i;
	}
	
	
    }
    

    return $best_ov;

}

sub mergeOverlappingSequences {
    my ($s1, $s2, $r) = @_;
    
    my @a = split //, $s2;

    if ($r == -1) {
	return $a[0] . $s1;
    } else {
	my $n = pop @a;
	return $s1 . $n;
    }

}


sub arraySum {
    my ($a_ref) = @_;

    my $n = scalar(@$a_ref);

    my $sum = 0;
    for (my $i=0; $i<$n; $i++) {
	$sum += $a_ref->[$i];
    }

    return $sum;

}


sub addArray2ToArray1 {
    my ($a_ref1, $a_ref2) = @_;

    my $n = @$a_ref1;

    for (my $i=0; $i<$n; $i++) {
	$a_ref1->[$i] = $a_ref1->[$i] + $a_ref2->[$i];	
    }
}


sub addArrays {
    my ($a_ref1, $a_ref2) = @_;
    
    my $a_ref3 = [];
    
    my $n = max(scalar(@$a_ref1), scalar(@$a_ref2));

    for (my $i=0; $i<$n; $i++) {
	$a_ref3->[$i] = $a_ref1->[$i] + $a_ref2->[$i];	
    }
    
    return $a_ref3;
}

#
#  calculate a distribution based on min, max, nbbins
#
sub getSpecifiedDistribution {
    
    my ($a_ref_corr, $i_nbbins, $i_min, $i_max) = @_;

    die "Cannot call Sets->getSpecificDistribution with a no sets of values\n"
	if (scalar(@$a_ref_corr) == 0);
    
    my @count = ();
    for (my $b=0; $b <$i_nbbins; $b++) {
	$count[$b] = 0;
    }

    # calculate the size of the bin from the 
    my $f_binsize = ($i_max - $i_min) / $i_nbbins;


    # fill each bin
    foreach my $c (@$a_ref_corr) {
	for (my $b=0; $b <$i_nbbins; $b++) {
	    if ( ($c > $i_min + $b * $f_binsize) && 
		 ($c < $i_min + ($b+1) * $f_binsize)) {
		$count[$b] ++;		
		last;
	    }
	}
    }
    
    my $total = scalar(@$a_ref_corr);

    # create an array of (bin middle, bin freq)
    my $t = 0;
    my @dist = ();
    for (my $b=0; $b <$i_nbbins; $b++) {	
	my @tmp =( $i_min + $b * $f_binsize + $f_binsize/2.0, $count[$b]/$total);
	
	push @dist, \@tmp;
    }

    return \@dist;    
}




sub isSubstringOf_2S {
    my ($s1, $s2, $diff) = @_;
	

    if (defined($diff)) {
	return 0 if ((length($s2)-length($s1)) != $diff);
    }
    
    my $ss1 = $s1;
    my $sc1 = Sets::getComplement($ss1);

    if (($s2 =~ /$ss1/) || ($s2 =~ /$sc1/)) {
	return 1;
    } else {
	return 0;
    }
    
}


sub isSubstringOf_1S {
    my ($s1, $s2, $diff) = @_;
	

    if (defined($diff)) {
	return 0 if ((length($s2)-length($s1)) != $diff);
    }
    
    my $ss1 = $s1;

    if ($s2 =~ /$ss1/) {
	return 1;
    } else {
	return 0;
    }
    
}





sub getDistribution {
    
    my ($a_ref_corr, $i_nbbins) = @_;

    die "Cannot call Sets->getDistribution with a no sets of values\n"
	if (scalar(@$a_ref_corr) == 0);
    my @count = ();
    for (my $b=0; $b <$i_nbbins; $b++) {
	$count[$b] = 0;
    }

    # calculate the size of the bin from the 
    my $f_binsize = 2.0 / $i_nbbins;

    foreach my $c (@$a_ref_corr) {
	for (my $b=0; $b <$i_nbbins; $b++) {
	    if ( ($c > -1 + $b * $f_binsize) && 
		 ($c < -1 + ($b+1) * $f_binsize)) {
		$count[$b] ++;
		
		last;
	    }
	}
    }
    
    my $total = scalar(@$a_ref_corr);


    
    # order the correlation set    
    
    my $t = 0;
    my @dist = ();

    for (my $b=0; $b <$i_nbbins; $b++) {
	
	my @tmp =( -1 + $b * $f_binsize + $f_binsize/2.0, $count[$b]/$total);
	
	push @dist, \@tmp;
    }

    return \@dist;    
}


#
# get the nb of identities between two strings
#
sub getNbMatches {
    my ($s1, $s2) = @_;
    
    my $l = length($s1);
    
    my @a = split //, $s1;
    my @b = split //, $s2;

    my $id = 0;
    for (my $i=0; $i<$l; $i++) {
	$id ++ if ($a[$i] eq $b[$i]);
    }

    #print "... comparing $s1 and $s2 ..\n";

    return $id;
}



#
# get the nb of identities between two strings
#
sub getNbMatches_mRNA {
    my ($s1, $s2, $rna) = @_;
    
    my $l = length($s1);
    
    my @a = split //, $s1;
    my @b = split //, $s2;

    my $id = 0;
    for (my $i=0; $i<$l; $i++) {
	$id ++ if ($a[$i] eq $b[$i]);

	if (
	    (($a[$i] eq 'A') && ($b[$i] eq 'G')) ||
	    (($a[$i] eq 'G') && ($b[$i] eq 'A')) ||
	    (($a[$i] eq 'T') && ($b[$i] eq 'C')) ||
	    (($a[$i] eq 'C') && ($b[$i] eq 'T')) ) {
	    
	    $$rna ++;
	}
    }

    

    #print "... comparing $s1 and $s2 ..\n";

    return $id;
}

#
#      
#   NNNNNNNNNNNNNNNNNNNNNNNNNNN
#   <--- decal ---> NNNNNNNNNNNNNNNNN
#
#
sub getBestOverlapDP { 
    my ($s1, $s2, $r_strand, $r_decal, $ma, $mm, $twostrand) = @_;

    my $l1 = length($s1);
    my $l2 = length($s2);
    my @a  = ();
    my $best_a = undef;
    my $i;
    my $j;
    my $s2c;
    
    
    
    my @a_s1 = split //, $s1;
    my @a_s2 = split //, $s2;
    
    
    #
    #  simply fills an array of added matches/mismatches
    #
    
    for ($i=0; $i<$l1; $i++) {
	$a[$i][0] = ($a_s1[$i] eq $a_s2[0]?$ma:$mm);
    }
    
    for ($j=0; $j<$l2; $j++) {
	$a[0][$j] = ($a_s1[0] eq $a_s2[$j]?$ma:$mm);
    }
    
    for ($i=1; $i<$l1; $i++) {
	for ($j=1; $j<$l2; $j++) {
	    $a[$i][$j] = $a[$i-1][$j-1] + ($a_s1[$i] eq $a_s2[$j]?$ma:$mm);
	}
    }


   # for ($i=1; $i<$l1; $i++) {
	#for ($j=1; $j<$l2; $j++) {
	 #   print "$a[$i][$j] ";
	#}
	#print "\n";
    #}
    
    #  at this point, look for the best diagonal
    
    $best_a = -100;
    
    
    for ($i=0; $i<$l1; $i++) {
	if ($a[$i][$l2 - 1] > $best_a) {
	    $best_a   = $a[$i][$l2 - 1];
	    $$r_decal   = - ($l2 - ($i + 1));
	    $$r_strand  = 1;
	}
    }
    
    for ($j=0; $j<$l2; $j++) {
	if ($a[$l1 - 1][$j] > $best_a) {
	    $best_a  = $a[$l1 - 1][$j];
	    $$r_decal  = $l1 - ($j + 1);
	    $$r_strand = 1;
	}
    }
    
    #return $best_a;
    

  if ($twostrand == 1) {
    
    # do the same, but reverse complement s2
    
    my $s2c = getComplement($s2);
    
    my @a_s2c = split //, $s2c;
    
    for ($i=0; $i<$l1; $i++) {
      $a[$i][0] = ($a_s1[$i] eq $a_s2c[0]?$ma:$mm);
    }
    
    for ($j=0; $j<$l2; $j++) {
      $a[0][$j] = ($a_s1[0] eq $a_s2c[$j]?$ma:$mm);
    }
    
    for ($i=1; $i<$l1; $i++) {
      
      for ($j=1; $j<$l2; $j++) {
	
	$a[$i][$j] = $a[$i-1][$j-1] + ($a_s1[$i] eq $a_s2c[$j]?$ma:$mm);
	
      }
      
    }
    
    #  at this point, look for the best diagonal
    
    
    for ($i=0; $i<$l1; $i++) {
      if ($a[$i][$l2 - 1] > $best_a) {
	$best_a   = $a[$i][$l2 - 1];
	$$r_decal   = - ($l2 - ($i + 1));
	$$r_strand  = 0;
      }
    }
    
    for ($j=0; $j<$l2; $j++) {
      if ($a[$l1 - 1][$j] > $best_a) {
	$best_a  = $a[$l1 - 1][$j];
	$$r_decal  = $l1 - ($j + 1);
	$$r_strand = 0;
      }
    }
    
    
  }

  return $best_a;
  
}


# get the best possible number of matches between two strings
sub getBestOverlap {
    my ($s1, $s2, $a_ref) = @_;
    
    #print "Comparing $s1 and $s2\n";
    
    
    # strand 1
    my $ov = 0;    
    for (my $i=0; $i<length($s1); $i++) {
	my $ss1 = substr($s1, $i); 
	my $ss2 = substr($s2, 0, length($s2) - $i);
	
	my $id  = getNbMatches($ss1, $ss2);
	if ($id > $ov) {
	    $ov = $id;
	    $a_ref->{STRAND} = 0;
	    $a_ref->{DECAL}  = $i;
	}
	
    } 

    for (my $i=0; $i<length($s2); $i++) {
	my $ss1 = substr($s2, $i); 
	my $ss2 = substr($s1, 0, length($s1) - $i);
	
	my $id  = getNbMatches($ss1, $ss2);
	if ($id > $ov) {
	    $ov = $id;
	    
	    $a_ref->{STRAND} = 0;
	    $a_ref->{DECAL}  = -$i;
	    
	}
	
    } 

   
   

    my $s1_cmp = getComplement($s1);
    #print "Comparing $s1_cmp and $s2\n";

    for (my $i=0; $i<length($s1_cmp); $i++) {
	my $ss1 = substr($s1_cmp, $i); 
	my $ss2 = substr($s2, 0, length($s2) - $i);
	my $id  = getNbMatches($ss1, $ss2);
	if ($id > $ov) {
	    $ov = $id;
	    
	    $a_ref->{STRAND} = 1;
	    $a_ref->{DECAL}  = -$i;
	}
	
    } 
    
     for (my $i=0; $i<length($s2); $i++) {
	my $ss1 = substr($s2, $i); 
	my $ss2 = substr($s1_cmp, 0, length($s1_cmp) - $i);
	my $id  = getNbMatches($ss1, $ss2);
	if ($id > $ov) {
	    $ov = $id;
	    
	    $a_ref->{STRAND} = 1;
	    $a_ref->{DECAL}  = $i;
	}
	
    } 
    
    return $ov;

}




# get the best possible number of matches between two strings

sub getBestOverlap_singlestrand {
    my ($s1, $s2, $a_ref) = @_;
    
    
    my $s1plus  = "X" x length($s2);
       $s1plus .= $s1;
       $s1plus .= "X" x length($s2);
    

    my $ov = 0;    
    for (my $i=0; $i<length($s1plus)-length($s2)+1; $i++) {

	my $s2plus  = "Y" x $i; 
	$s2plus    .= $s2;
	$s2plus    .= "Y" x (length($s1plus)-length($s2plus));

	#print "comparing $s1plus and $s2plus\n";

	my $rnamatches = 0;
	
	my $id  = getNbMatches_mRNA($s1plus, $s2plus, \$rnamatches);
	
	
#	print "RNA=$rnamatches\n";
	
	if ($id > $ov) {
	    $ov = $id;
	    $a_ref->{STRAND} = 0;
	    $a_ref->{DECAL}  = $i;
	    $a_ref->{MATCH}  = substr($s1plus, $i, length($s2));
	    $a_ref->{RNAMATCH}  = $rnamatches;
	    
	}
	
    } 

    

    return $ov;

}




sub sign {
  my ($m) = @_;

  if ($m == 0) {
    return 0;
  } elsif ($m > 0) {
    return 1;
  } else {
    return -1;
  }

}

sub getComplement {
    
    my ($str) = @_;
    
    my $l = length($str);
    
    my @s = split //, $str;

    my $c = "";
    for (my $i=$l-1; $i>=0; $i--) {
	
	my $d = "";
	
	if ($s[$i] eq 'A') {
	    $d = 'T';
	} elsif ($s[$i] eq 'T') {
	    $d = 'A';
	} elsif ($s[$i] eq 'G') {
	    $d = 'C';
	} elsif ($s[$i] eq 'C') {
	    $d = 'G';
	} elsif ($s[$i] eq '-') {
	    $d = '-';
	} elsif ($s[$i] eq '[') {
	    $d = ']';
	} elsif ($s[$i] eq ']') {
	    $d = '[';
	} elsif ($s[$i] eq '.') {
	    $d = '.';
	} else {
	    $d = 'N';
	}

	$c .= $d;
    }
    
    return $c;
   
}


sub getREComplement {
    
    my ($str) = @_;
    
    my $l = length($str);
    
    my @s = split //, $str;

    my $c = "";
    for (my $i=$l-1; $i>=0; $i--) {
	
	my $d = "";
	
	if ($s[$i] eq 'A') {
	    $d = 'T';
	} elsif ($s[$i] eq 'T') {
	    $d = 'A';
	} elsif ($s[$i] eq 'G') {
	    $d = 'C';
	} elsif ($s[$i] eq 'C') {
	    $d = 'G';
	} elsif ($s[$i] eq '[') {
	    $d = ']';
	} elsif ($s[$i] eq ']') {
	    $d = '[';
	} else {
	    $d = '.';
	}

	$c .= $d;
    }
    
    return $c;
   
}




sub getRNAComplement {
    
    my ($str) = @_;
    
    my $l = length($str);
    
    my @s = split //, $str;

    my $c = "";
    for (my $i=$l-1; $i>=0; $i--) {
	
	my $d = "";
	
	if ($s[$i] eq 'A') {
	    $d = 'U';
	} elsif ($s[$i] eq 'U') {
	    $d = 'A';
	} elsif ($s[$i] eq 'G') {
	    $d = 'C';
	} elsif ($s[$i] eq 'C') {
	    $d = 'G';
	} elsif ($s[$i] eq '-') {
	    $d = '-';
	} else {
	    $d = 'N';
	}

	$c .= $d;
    }
    
    return $c;
   
}


sub assemble_overlapping_fragments {

  my ($a_ref_segments) = @_;

  my %a = ();
  
  my $n = scalar( @$a_ref_segments );
  for (my $i=0; $i<$n; $i++) {
    for (my $j=0; $j<$n; $j++) {
      
      #next if ($i == $j);
      
      if (sequencesOverlap($a_ref_segments->[$i]->[0], $a_ref_segments->[$i]->[1], 
			  $a_ref_segments->[$j]->[0], $a_ref_segments->[$j]->[1])) {
	push @{ $a{ $i } }, $j; 
      }
    }
  }

  #foreach my $k (keys(%a)) {
  #  print "$k =>";
  #  foreach my $s (@{ $a{$k} }) {
  #    print " $s";
   # } 
   # print "\n";
  #}

  my $a_ref_c = get_connected_components( \%a );
  
  my @segs = ();

  foreach my $se (@$a_ref_c) {

   # print "Connected compo\n";
    
    my $i = shift @$se;
    my $s = $a_ref_segments->[$i]->[0];
    my $e = $a_ref_segments->[$i]->[1];

    foreach my $i (@$se) {
      if ($a_ref_segments->[$i]->[0] < $s) {
	$s = $a_ref_segments->[$i]->[0];
      }
      if ($a_ref_segments->[$i]->[1] > $e) {
	$e = $a_ref_segments->[$i]->[1];
      }
    }
    my @a_tmp = ($s, $e);
    push @segs, \@a_tmp;

    #print "s=$s e=$e\n";
  }

  return \@segs;

}



#
# returns 1 if sequences overlap, 0 otherwise
#   assumes that both sequences are oriented in the same sense
#   assumes BLAST sequence definition (1 .. n instead of 0..n)
sub sequencesOverlap {
    my ($s1, $e1, $s2, $e2) = @_;

    my $p1 = max( $s1,
		  $s2 );
    
    my $p2 = min( $e1,
		  $e2 );
    
    #my $li = $e1 - $s1 + 1;
    #my $lj = $e2 - $s2 + 1;
    my $ll = $p2 - $p1 + 1;
        
    if ( $ll > 0  ) {
	return $ll;
    } else {
	return 0;
    }
    
}

#
#  assume that the sequence overlap
#
sub getUnionForSequenceOverlap {
    my ($s1, $e1, $s2, $e2) = @_;

    my $p1 = min( $s1,
		  $s2 );
    
    my $p2 = max( $e1,
		  $e2 );
    
    return [$p1 , $p2];
        
    
}


sub getBondariesForSequenceOverlap {
    my ($s1, $e1, $s2, $e2) = @_;

    my $p1 = max( $s1,
		  $s2 );
    
    my $p2 = min( $e1,
		  $e2 );
    
    return [$p2 , $p1];
        
    
}



sub getSequencesOverlap {
    my ($s1, $e1, $s2, $e2) = @_;

    my $p1 = max( $s1,
		  $s2 );
    
    my $p2 = min( $e1,
		  $e2 );
    
    #my $li = $e1 - $s1 + 1;
    #my $lj = $e2 - $s2 + 1;
    my $ll = $p2 - $p1 + 1;
        
    return $ll;
    
}



#
# kmer generation
#
sub allkmers {
    
    my ($m) = @_;
        
    my @a = ();
    kmer(0, $m, "", \@a) ;
 
    return \@a;
}



sub kmer {
    
    my ($level, $maxlevel, $s, $r) = @_;
    
    if ($level == $maxlevel) {	
	push @$r, $s;
	return;
    }
    
    foreach my $n (("A", "C", "T", "G")) {
	kmer($level+1, $maxlevel, $s . $n, $r);
    }
}


#
# remove complemets : takes a list of kmers, output the same list with complements removed
#
sub removeComplements {
    my ($a_ref) = @_;
    
    my @a = ();
    
    # store each kmer into a hash
    my %index = ();
    foreach my $s1 (@$a_ref) {
	my $s2 = getComplement($s1);
	
	# if none of them is there
	if (!$index{$s1} && !$index{$s2}) {
	    $index{$s1} = 1;
	    
	    push @a, $s1;
	}
	
	
    }

    return \@a;
}

#
#  count the nucleotides A, C, G, T
#
sub countNucleotides {

    my ($s) = @_;
    
    my %H = ();
    my @nt = split //, $s;
    foreach my $bp (@nt) {
	$H{ $bp } ++;
    }

    return \%H;
}

sub removeLowComplexityKmers {
    my ($a_ref, $n) = @_;
    
    my @a = ();
    
    foreach my $s (@$a_ref) {

	my @nt = split //, $s;
	my $na = 0;
	my $nt = 0;
	foreach my $bp (@nt) {
	    $na ++ if ($bp eq 'A');
	    $nt ++ if ($bp eq 'T');
	}
	
	next if ($na >= $n);
	next if ($nt >= $n);
	
	push @a, $s;
    }
    
    return \@a;
    
}

# returns the number of ATs
sub getPolyAT {
    my ($s) = @_;

    
    my @nt = split //, $s;
    my $na = 0;
    my $nt = 0;
    foreach my $bp (@nt) {
	$na ++ if ($bp eq 'A');
	$nt ++ if ($bp eq 'T');
    }

    return $na+$nt;
	
    
}


# returns the number of ATs
sub getNbLetterOccurrences {
    my ($s, $l) = @_;

    
    my @nt = split //, $s;
    my $nt = 0;
    foreach my $bp (@nt) {
	$nt ++ if ($bp eq $l);
    }

    return $nt;
	
    
}


sub dirname {
  my ($s) = @_;
  my ($p) = $s =~ /(.+)\/[^\/]+$/;
  if (!defined($p)) {
    $p = ".";
  } 
  return $p;
}


sub basename {
    my ($s) = @_;
    
    my ($p) = $s =~ /\/{0,1}?([^\/]+?)\.[^\.]+$/;

    return $p;
}


sub filename {
    my ($s) = @_;
    
    my ($p) = $s =~ /\/{0,1}?([^\/]+)$/;
    #my ($p) = $s =~ /([^\/]+)$/;

    return $p;
}


#
# remove the .1 in the WORM orf name
#
sub wormRenameSet {    
    my ($a_ref) = @_;
    my @a_tmp = ();    
    foreach my $o (@$a_ref) {
	my $oo = $o;
	if ($oo =~ /\.\d+\.\d+/) {
	    $oo =~ s/\.1$//;
	}	
	push @a_tmp, $oo;
    }    
    return \@a_tmp;
}


#
# print a set
#
sub printSet {
    my ($a_ref) = @_;
    foreach my $c (@$a_ref) {
	print "$c\n";
    }

}


sub printSetSep {
    my ($a_ref, $s) = @_;
    my @a  =();
    foreach my $c (@$a_ref) {
	push @a, $c;
    }
    print join($s, @a);
    print "\n";

}


#
# print a hash
#
sub printHash {
    my ($h, $s) = @_;
    
    if (!defined($s)) {
	$s = " => ";
    }
    
    foreach my $k (keys(%$h)) {
	print "$k$s" . $h->{$k} . "\n";
    }
    
}



sub writeSet {
    my ($a_ref, $file) = @_;

    open OO, ">$file" or die "Sets::writeSet : cannot open file \"$file\" ..\n";
    
    foreach my $c (@$a_ref) {

	print OO "$c\n";
	
    }

    close OO;
    
    
}


sub writeTable{
    my ($a_ref, $file) = @_;

    open OO, ">$file" or die "Sets::writeSet : cannot open file \"$file\" ..\n";
    
    foreach my $c (@$a_ref) {

	print OO join("\t", @$c); print OO "\n";
	
    }

    close OO;
}

sub countSymbols {
    my ($a_ref) = @_;

    my %HASH = ();
    foreach my $r (@$a_ref) {
	
        $HASH { $r } ++;
	
    }
    
    my $a_ref_sorted = Sets::getArrayKVFromHash(\%HASH, 1);

    return $a_ref_sorted;
}


sub writeSetToFile {
    my ($a_ref, $file) = @_;

    open OO, ">$file" or die "Sets::writeSet : cannot open file \"$file\" ..\n";
    
    foreach my $c (@$a_ref) {

	print OO "$c\n";
	
    }

    close OO;
    
    
}


sub writeDistribution {
    my ($a_ref, $file) = @_;

    open OO, ">$file";
    
    foreach my $c (@$a_ref) {

	print OO "$c->[0]\t$c->[1]\n";
	
    }

    close OO;
    
    
}


#
#  
#
sub in_array {
    my $val = shift(@_);
 
    foreach my $elem (@_) {
        if($val eq $elem) {
            return 1;
        }
    }

    return 0;
}


sub getIntersectionSet {
    my ($a_ref1, $a_ref2) = @_;
    return getOverlapSet($a_ref1, $a_ref2);
}

sub getIntersection {
    my ($a_ref1, $a_ref2) = @_;
    return getOverlapSet($a_ref1, $a_ref2);
}

#
# get the overlap set between two sets
#
sub getOverlapSet {
    my ($a_ref1, $a_ref2) = @_;
    

    my @a_overlap = ();
    my %h_there   = ();
    foreach my $k (@$a_ref1) {
	$h_there{$k} = 1;
    }

    foreach my $k (@$a_ref2) {
	push @a_overlap, $k if ($h_there{$k} == 1);
    }

    return \@a_overlap;
}



#
# get the overlap set between two sets
#
sub getNotIntersection {
    my ($a_ref1, $a_ref2, $n) = @_;

    my @a_otherset = ();
    my %h_there1   = ();
    my %h_there2   = ();

    #print "fuck?\n";

    #print scalar(@$a_ref2);
    
    foreach my $k (@$a_ref1) {
	$h_there1{$k} = 1;
    }

    foreach my $k (@$a_ref2) {
	$h_there2{$k} = 1;
    }

    if ($n == 0) {

	# return elements in 1 but not in 2
	foreach my $k (@$a_ref1) {
	    push @a_otherset, $k if (!defined($h_there2{$k}));
	}
    } else {
	foreach my $k (@$a_ref2) {
	    push @a_otherset, $k if (!defined($h_there1{$k}));
	}
    }

    return \@a_otherset;
}




sub getOverlapSize {

    my ($a_ref1, $a_ref2) = @_;
    
    my $a_ref3  = getOverlapSet($a_ref1, $a_ref2);

    return scalar(@$a_ref3);
}

#
#  get the union between two sets
#
sub getUnionSet {
    
    my ($a_ref1, $a_ref2) = @_;
    

 
    my %h_there   = ();

    foreach my $k (@$a_ref1) {
	$h_there{$k} = 1;
    }

    foreach my $k (@$a_ref2) {
	$h_there{$k} = 1;
    }

    my @a_union = keys( %h_there );
    
    return \@a_union;

}

#
#  transform a set into a SQL set
#
sub SetToSQLSet {
    
    my ($a_ref, $delim) = @_;
    
    my @a_tmp = @$a_ref;

    foreach my $r (@a_tmp) {
	$r = "$delim$r$delim";
    }
    
    my  $s = join (",", @a_tmp);

    return "($s)";
    
}

#
#  SQL ref set to set 
#
sub SQLRefToSet {
    my ($a_ref, $k) = @_;
    
    my @a_orfs = ();
    foreach my $d_ref (@$a_ref) {
	push @a_orfs, $d_ref->{$k};
    }
    
    return \@a_orfs;
}

#
#  SQL ref set to simple hash
#
sub SQLRefToSimpleHash {
    my ($a_ref, $k, $v) = @_;
    
    my %h = ();
    foreach my $d_ref (@$a_ref) {
	$h{$d_ref->{$k}} = $d_ref->{$v};
    }
    
    return \%h;
}


#
#  SQL ref set to simple hash
#
sub SQLRefToIndex {
    my ($a_ref, $k) = @_;
    
    my %h = ();
    foreach my $d_ref (@$a_ref) {
	$h{$d_ref->{$k}} = $d_ref;
    }
    
    return \%h;
}




# read a set from disk
# $d = 1 => remove duplicates
sub readSet {
    my ($file, $d) = @_;
    
    open IN, $file or die "$file cannot be opened\n";
    my @a_tmp = ();
    while (my $l = <IN>) {
	chomp $l;
	push @a_tmp, $l if ($l ne "");
    }
    close IN;

    if (defined($d) && ($d == 1)) {
	my $r = removeDuplicates(\@a_tmp);
	return $r;
    } else {
	return \@a_tmp;
    }
}


#
# get the size of a set (nb of lines)
#
sub getSetSize {
    my ($file)  = @_;

    my $a_ref = readSet($file);
    
    return scalar(@$a_ref);
}


#
# get an index from a set in a file
#
sub getIndex {
    my ($f) = @_;

    my %i = ();
    
    open IN, $f;
    while (my $l = <IN>) {
	chomp $l;

	$i{$l} = 1;
    }
    close IN;

    return \%i;
}


#
# get an index from a set in a file
#
sub getIndexFromTableColumn {
    my ($f, $p) = @_;

    my %i = ();
    
    open IN, $f;
    while (my $l = <IN>) {
	chomp $l;

	my @a = split /\t/, $l, -1;
	
	$i{$a[$p]} = 1;
    }
    close IN;

    return \%i;
}


#
# get an index from a set in a file
#

sub getIndexKV {
    my ($f, $x, $y) = @_;

    my %i = ();
    
    open IN, $f;
    while (my $l = <IN>) {
	chomp $l;
	my @a = split /\t/, $l, -1;
	$i{$a[$x]} = $a[$y];
    }
    close IN;

    return \%i;
}


#
# get an index from a set in a file
#
sub getIncreasingIndex {
    my ($f) = @_;

    my %i = ();
    
    open IN, $f; my $i = 1;
    while (my $l = <IN>) {
	chomp $l;

	$i{$l} = $i;

	$i++;
    }
    close IN;

    return \%i;
}



#
# get an index from a set ref
#
sub  getIndexFromArrayRef {
    my ($a_ref) = @_;

    my %i = ();
    
    foreach my $r (@$a_ref) {

	$i{$r} = 1;
    }



    return \%i;
}


sub  getIncreasingIndexFromArrayRef {
    my ($a_ref) = @_;

    my %i = ();
    my $u = 0;
    foreach my $r (@$a_ref) {
	$i{$r} = $u;
	$u++;
    }
    return \%i;
}



#
# get an index from a set ref
#
sub  getIndexRefFromArrayRef {
    my ($a_ref, $ab) = @_;

    my %i = ();
    
    foreach my $r (@$a_ref) {

	$i{ $r->[$ab] } = $r;
    }



    return \%i;
}



#
#  get a bi-array from a hash
#
sub getArrayKVFromHash {
    my ($h_ref, $sorted) = @_;

    my @ar = ();
    foreach my $k (keys(%$h_ref)) {
	
	my @a_tmp = ($k, $h_ref->{$k});
	
	push @ar, \@a_tmp;
    }

    if ($sorted == 1) {
	@ar = sort { $a->[1] <=> $b->[1] } @ar;
    }

    return \@ar;


}



sub getMultipleIndexKV {
    my ($a_ref, $k, $v) = @_;

    my $h_ref_regulatory_sites = undef;
    foreach my $r (@$a_ref) {
	push @{ $h_ref_regulatory_sites->{$r->[$k]} }, $r->[$v];
    }

    return $h_ref_regulatory_sites;
    

}


sub removeDuplicates {
    my ($a_ref) = @_;
    
    my @b = ();
    
    # inverted inex
    my %ix = ();
    foreach my $aa (@$a_ref) {
	push @b, $aa if (!$ix{$aa});
	$ix{$aa} = 1;
    }

    #print join("\n", @b) . "\n";
    
    return \@b;
    
}

#
#  return positions starting at 0 ?
#
sub getREMotifPositions {
    my ($s_motif, $s_sequence) = @_;
    my @a_positions =();
    while ($s_sequence =~/$s_motif/ig){
        push (@a_positions, pos($s_sequence) - length($&));
    }

    my $mc = getComplement($s_motif);
    
    while ($s_sequence =~/$mc/ig){
        push (@a_positions, pos($s_sequence) - length($&));
    }
        
    return removeDuplicates(\@a_positions);
}


#
#  same as above, but also returns orientation
#
sub getREMotifPositionsOrientations {

  my ($s_motif, $s_sequence) = @_;
  
  my %h_pos = ();

  my @a_positions =();
  while ($s_sequence =~/$s_motif/ig){
    my $p = pos($s_sequence) - length($&);
    my @t = ($p, 1); 
    push @a_positions, \@t;
  }
  
  my $mc = getComplement($s_motif);    
  while ($s_sequence =~/$mc/ig){
    my $p = pos($s_sequence) - length($&);
    my @t = ($p, -1); 
    push @a_positions, \@t;
  }
  
  return \@a_positions;

}





sub getREMotifPositions_singlestrand {
    my ($s_motif, $s_sequence) = @_;
    my @a_positions =();
    while ($s_sequence =~/$s_motif/ig){
        push (@a_positions, pos($s_sequence) - length($&));
    }

        
    return \@a_positions;
}



sub getRandomElement {
    my ($a_ref) = @_;

    my $i_size   = scalar(@$a_ref);  
    my $i = int(rand($i_size));

    return $a_ref->[$i];
}


#
# randomly sample $n elements from a set (with resampling)
# 
sub sampleSet {
    
    my ($n, $a_ref) = @_;

    srand;
    
    my $i_size   = scalar(@$a_ref);    
    
    my @in       = ();

    my @a_tmp = ();
    
    for (my $a=0; $a<$n; $a++) {

        my $i = int(rand($i_size));
    
	#if ($r && (in_array($i,

        my $o = $a_ref->[$i];
            
	push @a_tmp, $a_ref->[$i];
    
    }

    return \@a_tmp;
}


sub sampleSetNoRep {
    
    my ($n, $a_ref) = @_;


    
    my $i_size   = scalar(@$a_ref);    
    
    my @in       = ();

    my @a_tmp = ();
    
    for (my $a=0; $a<$n; $a++) {

        my $i = int(rand($i_size));
	while (Sets::in_array($a_ref->[$i], @a_tmp)) {
	    $i = int(rand($i_size));
	} 

	#if ($r && (in_array($i,

        my $o = $a_ref->[$i];
            
	push @a_tmp, $a_ref->[$i];
    
    }

    return \@a_tmp;
}


#
# returns a temptp file name
#
sub getTempFile {
    my ($prefix, $suffix) = @_;
    my $d = time;
    my $h = $prefix . int(rand(1000)) . shuffle_seq($d);
    if (defined($suffix)) {
      $h .= $suffix;
    }
    return $h;
}


sub saveToTempFile {
    my ($s) = @_;

    my $tmp = getTempFile("/tmp/tchutchu");
    
    open HOUTCH, ">$tmp";
    print HOUTCH $s;
    close HOUTCH;
    
    return $tmp;
}

sub getRandomString {
    my ($prefix) = @_;
	
   
    my $d = time;
    
    my $h = $prefix . "." . int(rand(1000)) . "." . shuffle_seq($d);
     
    #print "$h\n";
    return $h;
    
}



sub log2 {

    my ($l) = @_;
    
    return log($l) / log(2);
}




sub log10 {

    my ($l) = @_;
    
    return log($l) / log(10);
}




sub shuffle_array {
    my ($a_old) = @_;
    
    my @a_new = ();
    
    for( @$a_old ){
        my $r = rand @a_new+1;
        push(@a_new,$a_new[$r]);
        $a_new[$r] = $_;
    }

    return \@a_new;
}


#
#  shuffle a sequence
#
sub shuffle_seq {
    my ($s) = @_;
    
    my @a = split //, $s;    
    my $a_ref = shuffle_array(\@a);
    return join("", @$a_ref);
}


#
#  shuffle k-mers
#
sub shuffle_seq_kmers {
    
    my ($s, $k) = @_;
    
    my @a = ();
    my $p = undef;
    my $pat = "(.{$k})";
    while ($s =~ /$pat/gi) {
	#print "$1\n";
	push @a, $1;

	$p =  pos($s); 
    }

    
    push @a, substr($s, $p) if ($p != length($s));
    
    #print join("\n", @a); print "\n";

    
    my $a_ref = shuffle_array(\@a);
    return join("", @$a_ref);

}


sub readKmersIdx {
    my ($f) = @_;
    my %idx = ();
    open INICH, $f;
    while (my $l = <INICH>) {
        #print "$l";
        chomp $l;
        my @a = split /\t/, $l;
        $idx{$a[0]} = \@a;
    }
    close INICH;
    
    return \%idx;
}    


sub readKmers {
    my ($f) = @_;
    my @array = ();
    open INICH, $f or die "cannot open $f\n";
    while (my $l = <INICH>) {
	chomp $l;
	my @a = split /\t/, $l;
        push @array, \@a;
    }
    close INICH;
    
    return \@array;
}    



sub getArray {
    my ($f) = @_;
    my @array = ();
    open INICH, $f or die "cannot open $f\n";
    while (my $l = <INICH>) {
	chomp $l;
	my @a = split /\t/, $l, -1;
        push @array, \@a;
    }
    close INICH;
    
    return \@array;
}    



sub max {
    my ($a, $b) = @_;
    return ($a>$b?$a:$b);
}



sub minInArray {
    my ($r) = @_;

    my @aa = sort { $a <=> $b } @$r;

    
    return $aa[ 0 ];
}


sub maxInArray {
    my ($r) = @_;

    my @aa = sort { $a <=> $b } @$r;

    
    return $aa[ $#aa ];
}


sub indexMaxInArray {
  my ($r) = @_;
  my $aa = order($r, 0);
  return $aa->[0];
}


sub percentile {

    my ($r, $p) = @_;
    my @aa = sort { $a <=> $b } @$r;

    my $n  = scalar(@aa);

    my $idx = $p * $n;

    my $idx_do = int($idx + 0.5)    ;
    if ($idx_do > $n-1) {
	$idx_do =  $n-1;
    }
    #my $idx_up = int($idx) + 1;

    
    return $aa[ $idx_do ]; # + $aa[ $idx_up ]) / 2;
}




sub min {
    my ($a, $b) = @_;
    return ($a<$b?$a:$b);
}

#
#  generate a random symbol from a HASH
#
sub generateRandomSymbol {
    my ($h_ref) = @_;
    
    my @a_cum = ();
    my $i     = 1;
    my @a     = keys(%$h_ref);
    
    $a_cum[0] = 0.0;
    foreach my $l1 (@a) {
	$a_cum[$i] = $a_cum[$i - 1] + $h_ref->{ $l1 };
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
    
    return $n;


}



#
#  read an Ensembl ortholog prediction
#    returns an index GENENAME -> ARRAY
#    only keeps the best ortholgs !
sub readEnsemblOrthologs {
    
    my ($file) = @_;
    
    my %DATA = ();

    open IN1, $file or die "cannot read ortholog file ..\n";
    while (my $l = <IN1>) {
	chomp $l;
	next if ($l =~ /Chromosome/);

	my @a_tmp = split /\t/, $l;

	# next if ortholog already there and score does not improve
	next if (defined($DATA{$a_tmp[4]} && ($DATA{$a_tmp[4]}->[9] > $a_tmp[9]))); 

	$DATA{$a_tmp[4]} = \@a_tmp;
    }

    close IN1;

    return \%DATA;

}

#
#  
#
sub getReciprocalOrthologs {
    my ($h_ref1, $h_ref2) = @_;

    my @a = ();
    
    foreach my $k1 (keys(%$h_ref1)) {

	
	# get the best ortholog in species 2
	my $bo2 = $h_ref1->{$k1}->[5];
	
	# get the best ortholog in species 1
	my $bo1 = $h_ref2->{$bo2}->[5];
	
	if ($k1 eq $bo1) {
	    my @a_tmp  = ($k1, $bo2);
	    push @a, \@a_tmp;
	}
	
    }
    
    return \@a;
}

#
#  calculate the standard deviation of a series of values
#
sub stddev {
    my ($a_ref) = @_;
    
    my $m = scalar(@$a_ref);

    return -1 if ($m <= 1);

    # calculate the sums
    my $sum   = 0.0;
    my $sum_2 = 0.0;
    for (my $i=0; $i<$m; $i++) {
	$sum   += $a_ref->[$i];
	$sum_2 += $a_ref->[$i] * $a_ref->[$i];
    }
   
    # calculate the standard deviation
    my $std = sqrt( ($sum_2 - $sum * $sum / $m ) / ( $m - 1 )); 
  
    return $std;

}


sub average {
    my ($a_ref) = @_;
    
    my $m = scalar(@$a_ref);

    # calculate the sums
    my $sum   = 0.0;
    
    for (my $i=0; $i<$m; $i++) {
	$sum   += $a_ref->[$i];
	
    }
   
    
  
    return $sum / $m;
    
    
}



# compute the median value of an array
sub median {
    
    my ($a_ref) = @_;
    
    my $n     = scalar(@$a_ref);
    my @a_tmp = sort {$a <=> $b} @$a_ref;
    
    my $nm    = int($n / 2)-1;  # floor - 1
    
    if ($n % 2 == 1)  {
        return $a_tmp[$nm + 1];
    } else {
        return ($a_tmp[$nm] + $a_tmp[$nm + 1]) / 2;
    }

    
}


sub pearson {
    
    my ($a_ref1, $a_ref2) = @_;

    my $m1 = scalar(@$a_ref1);
    my $m2 = scalar(@$a_ref2);
    
    # calculate the average
    my $sum1   = 0.0;
    my $sum2   = 0.0;

    my $actual_m = 0;
    
    for (my $i=0; $i<$m1; $i++) {

	next if (($a_ref1->[$i] eq "") || ($a_ref2->[$i] eq ""));

	#print $a_ref1->[$i] . " " . $a_ref2->[$i] . "\n";

        $sum1   += $a_ref1->[$i];
        $sum2   += $a_ref2->[$i];

	$actual_m++;
			  
    }

    

    my $avg1 = $sum1 / $actual_m;
    my $avg2 = $sum2 / $actual_m;

    # calc the Pearson correlation
    
    my $a = 0.0;
    my $b = 0.0;
    my $c = 0.0;
    for (my $i=0; $i<$m1; $i++) {

	next if (($a_ref1->[$i] eq "") || ($a_ref2->[$i] eq ""));


   
	$a += ($a_ref1->[$i] - $avg1) * ($a_ref2->[$i] - $avg2);
        $b += ($a_ref1->[$i] - $avg1) * ($a_ref1->[$i] - $avg1);
        $c += ($a_ref2->[$i] - $avg2) * ($a_ref2->[$i] - $avg2);

    }
    
    if (($b == 0.0) || ($c == 0.0)) {
	return 0.0;
    }
    
    return $a / sqrt ( $b * $c );
}


#
#   get the euclidean distance
#
sub euclidean {
    
    my ($a_ref1, $a_ref2) = @_;

    my $m1 = scalar(@$a_ref1);
    my $m2 = scalar(@$a_ref2);
    
    die "Both vectors must have the same size .\n" if ($m1 != $m2);

    my $sum    = 0.0;
    
    for (my $i=0; $i<$m1; $i++) {
        $sum   += ( $a_ref1->[$i] - $a_ref2->[$i] ) * ( $a_ref1->[$i] - $a_ref2->[$i] );
    }
    
    return sqrt($sum);
}


sub jacquard {
    my ($r1, $r2) = @_;

    my $a_int = getOverlapSet($r1, $r2);
    my $a_uni = getUnionSet($r1, $r2);

    my $ja   = scalar(@$a_int) / scalar(@$a_uni);

    return $ja;
}


sub Sets::size {
    my ($r1) = @_;

    return scalar(@$r1);
}

sub getFiles {
    my ($ls) = @_;



    my $s_list = `ls $ls`;

    #print "$s_list";

    my @a = split /\n/, $s_list;

    return \@a;


}





sub getManyFiles {
    my ($ls) = @_;
    my $s_list = `find /home/olly/DATA/YEASTS/ALIGNMENTS/ -name "*.aln"`;
    my @a = split /\n/, $s_list;
    return \@a;
}

sub file2txt {
    my ($f) = @_;

    open INF, $f or die "cannot open $f ..\n";
    my @a = <INF>;
    close INF;
    
    return join("", @a);
}


sub getNumberDifferencesBetweenKmers {
    my ($s1, $s2) = @_;

    my $l = length($s1);
    my @a = split //, $s1;
    my @b = split //, $s2;
    my @c = split //, getComplement($s2);
    
    #  first sense
    my $cnt1 = 0;
    for (my $k=0; $k<$l; $k++) {
	if ($a[$k] ne $b[$k]) {
	    $cnt1 ++;
	}
    }

    #  second sense
    my $cnt2 = 0;
    for (my $k=0; $k<$l; $k++) {
	if ($a[$k] ne $c[$k]) {
	    $cnt2 ++;
	}
    }
    
    return min($cnt1, $cnt2);
}




sub getSequencesIdentity {
    my ($s1, $s2) = @_;

    my $l = length($s1);
    my @a = split //, $s1;
    my @b = split //, $s2;
    
    #  first sense
    my $cnt1 = 0;
    for (my $k=0; $k<$l; $k++) {
	if ($a[$k] eq $b[$k]) {
	    $cnt1 ++;
	}
    }
    
    return $cnt1/$l;
}




sub trim {
    my ($s) = @_;

    my $ss = $s;

    $ss =~ s/[\ \t]//g; 
    
    return $ss;
}


sub translate {
    my ($s) = @_;
    
    my %CODON_TABLE = (
   TCA => 'S',TCG => 'S',TCC => 'S',TCT => 'S',
   TTT => 'F',TTC => 'F',TTA => 'L',TTG => 'L',
   TAT => 'Y',TAC => 'Y',TAA => '*',TAG => '*',
   TGT => 'C',TGC => 'C',TGA => '*',TGG => 'W',
   CTA => 'L',CTG => 'L',CTC => 'L',CTT => 'L',
   CCA => 'P',CCG => 'P',CCC => 'P',CCT => 'P',
   CAT => 'H',CAC => 'H',CAA => 'Q',CAG => 'Q',
   CGA => 'R',CGG => 'R',CGC => 'R',CGT => 'R',
   ATT => 'I',ATC => 'I',ATA => 'I',ATG => 'M',
   ACA => 'T',ACG => 'T',ACC => 'T',ACT => 'T',
   AAT => 'N',AAC => 'N',AAA => 'K',AAG => 'K',
   AGT => 'S',AGC => 'S',AGA => 'R',AGG => 'R',
   GTA => 'V',GTG => 'V',GTC => 'V',GTT => 'V',
   GCA => 'A',GCG => 'A',GCC => 'A',GCT => 'A',
   GAT => 'D',GAC => 'D',GAA => 'E',GAG => 'E',
   GGA => 'G',GGG => 'G',GGC => 'G',GGT => 'G');
    
    my @a = split //, $s;
    my $l = length($s);
    
    my $t = "";
    for (my $i=0; $i<$l; $i+=3) {
	my $c = $CODON_TABLE { $a[$i] . $a[$i+1]  . $a[$i+2] };
	$t .= ($c?$c:"?");
    }

    return $t;
}


sub getArrayOfCodons {
    my ($s) = @_;

    my @a = split //, $s;
    my $l = length($s);
    
    my @cs = ();
    for (my $i=0; $i<$l; $i+=3) {
	my $c =  $a[$i] . $a[$i+1]  . $a[$i+2];
	push @cs, $c;
    }
    
    return \@cs;

}


sub maskExons {
    my ($s, $a_ref_exons_boundaries, $c) = @_;

    my $seq = $s;
    foreach my $e (@$a_ref_exons_boundaries) {
	my $str = substr($seq, $e->[0], $e->[1]-$e->[0]+1);
	substr($seq, $e->[0], $e->[1]-$e->[0]+1) = $c x length($str); 
    }
    return $seq;
}


sub execRecompare {
    my ($re, $f1, $f2, $nb, $ts) = @_;

    my $todo = "/home/olly/PROGRAMS/FASTCOMPARE/recompare -fasta1 $f1 -fasta2 $f2 -re \"$re\" -nbgenes $nb -twostrand $ts -out toto.txt";
    
    #print $todo;
    
    my $txt  = `$todo`;  
    
    my ($score) = $txt =~ /\= ([\-\d\.]+)$/; 
    my ($ov)    = $txt =~ /overlap\=(\d+?),/;

    #print "sc=$score ov=$ov\n";

    my @a = ($score, $ov);

    return \@a;

}


sub getMaxCountFromHash {
  my ($h_ref) = @_;

  my $v_max  = undef;
  my $k_best = undef;
  while (my ($k, $v) = each(%$h_ref)) {
    
    if (!defined($k_best)) {
      $v_max  = $v;
      $k_best = $k;
    }

    elsif ($v > $v_max) {
      $v_max  = $v;
      $k_best = $k;
    }
    
  }

  return $k_best;
}

sub align_nt_on_aa_sequence {

  my ($s_aa1_aln, $s_nt1) = @_;
  
  my $s_nt1_aln = "";

  print "AA=$s_aa1_aln, NT=$s_nt1\n";

  my $c = $s_aa1_aln; $c =~ s/\-//g;
  print length($c); print "\n";
  print length($s_nt1); print "\n";

  # align the nt sequence to it
  my @a = split //, $s_aa1_aln;
  my $i = 0;
    

  foreach my $aa (@a) {
    
    #print " ** $aa ** \n";
    
    if ($aa eq "-") {
      $s_nt1_aln .= "---";
      print "$aa => ---\n";      
    } else {
      my $codon = substr($s_nt1, $i, 3);
      print "$aa => ". $codon . " (" . translate($codon) . ")\n";
      $s_nt1_aln .= substr($s_nt1, $i, 3);
      $i += 3;
    }
  }
 
  my $left = substr($s_nt1, $i);

  die "pb, left = $left\n" if ($left ne "");
  
  return $s_nt1_aln;

}


1;
