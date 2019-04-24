package Matrix2EPS;

=head1 NAME

 Matrix2EPS

=head1 SYNOPSIS

 use strict;
 use Matrix2EPS;

 my $ma = Matrix2EPS->new;
 $ma->setOutputFileName("out.eps");
 $ma->setMin(-3);
 $ma->setMax(3);
 $ma->setH(10); 
 $ma->setW(6);

 $ma->setColMap("HEATMAPS/cmap2.txt");
 $ma->loadMatrix($ARGV[0], 0, 0);
 $ma->loadIclustPartition($ARGV[1]);
 $ma->reorderUsingPartition(1,1);

 $ma->setMatrix(\@M, \@motifs, \@regions);
 $ma->clusterRows();
 $ma->clusterColumns("euclidean");

 $ma->setXbase(80);
 $ma->setYbase(180);

 $ma->setFont("Arial");
 $ma->setRownameFontSize(10);
 $ma->setHeaderFontSize(6);
 $ma->colnameAngle(45);

 $ma->loadColumnsPartition($ARGV[1]);


 $ma->drawBoxUp($val); # draw frame


 $ma->draw();
 $ma->addColumnClusters();
 $ma->addClusters();
 $ma->drawScale(20,180);
 
 $ma->output();
 $ma->pdfify();

=cut


use lib "$ENV{HEATMAPDIR}/PostScript-Simple-0.07/lib";
use lib "$ENV{HEATMAPDIR}";
use Table;
use Sets;
use PostScript::Simple;
use strict;

sub new {
  my ($class)       = @_;

  my $self        = {};

  $self->{MATRIX}   = undef;
  $self->{HEADER}   = undef;
  $self->{ROWNAMES} = undef;
  $self->{OUTFILE}  = "out.eps";
  $self->{COLMAP}   = undef;
  $self->{XBASE}    = 55;
  $self->{YBASE}    = 55;
  $self->{H}        = 5;
  $self->{W}        = 5;
  $self->{P}        = undef;
  $self->{VERBOSE}  = 1;
  $self->{PARTITION}= undef;
  $self->{MIN}               = -3;
  $self->{MAX}               =  3;
  $self->{MINTEXT}           = "down-regulation";
  $self->{MAXTEXT}           = "up-regulation";  
  $self->{DRAWFRAMES}        = 0;
  $self->{RFONTSIZE}         = 5;
  $self->{HFONTSIZE}         = 5;
  $self->{COLUMNS_PARTITION} = undef;
  $self->{DRAWROWNAMES}      = 'left';
  $self->{XRIGHT}            = 50;
  $self->{NEG}               = 0;
  $self->{COLDESC}           = undef;
  $self->{ROWDESC}           = undef;
  $self->{FONT}              = "Times";
  $self->{ARROWS}            = undef;
  $self->{DRAWHEADER}        = 1;
  $self->{ALGOCLUST}         = "avg";
  $self->{COLDRAWMOTIFS}     = 0;
  $self->{ROWDRAWMOTIFS}     = 0;
  $self->{WEBLOGODIR}        = "$ENV{HEATMAPDIR}/";
  $self->{DRAWSIGNIFBOXES}   = 0;
  $self->{DRAWBOX_UP_T}      = undef;
  $self->{DRAWBOX_DO_T}      = undef;
  $self->{COLNAME_ANGLE}     = 90;
  $self->{DRAWP}             = 1;
  $self->{MINABSLP}          = 2;
  $self->{DRAWHORIZSCALE}    = 0;
  $self->{ONLYSIGNIF}        = 0;
  $self->{ONLYPOS}           = 0;
  $self->{LPINCREMENTS}      = 1;
  $self->{SCALERES}          = 50;
  $self->{SCALEFONTSIZE}     = 14;
  $self->{RIGHTTABLE}        = undef;
  $self->{COLTOCHANGESIGN}   = undef;

  bless  $self;
  return $self;
}


sub changeNamedColSign {
  my ($self, $f) = @_;
  $self->{COLTOCHANGESIGN} = $f;
}

sub setScaleFontSize {
  my ($self, $f) = @_;
  $self->{SCALEFONTSIZE} = $f;
}

sub setScaleRes {
  my ($self, $f) = @_;
  $self->{SCALERES} = $f;
}


sub setMinAbsLP {
  my ($self, $f) = @_;
  $self->{MINABSLP} = $f;
}

sub setLPincrements {
  my ($self, $f) = @_;
  $self->{LPINCREMENTS} = $f;
} 

sub setFont {
  my ($self, $f) = @_;
  $self->{FONT} = $f;
}

sub colnameAngle {
  my ($self, $f) = @_;
  $self->{COLNAME_ANGLE} = $f;
}

sub drawSignifBoxes {
  my ($self, $f) = @_;
  $self->{DRAWSIGNIFBOXES} = $f;    
}

sub drawOnlySignif {
  my ($self, $f) = @_;
  $self->{ONLYSIGNIF} = $f;
}

sub drawOnlyPos {
  my ($self, $f) = @_;
  $self->{ONLYPOS} = $f;
}


sub drawBoxUp {
  my ($self, $f) = @_;
  $self->{DRAWBOX_UP_T} = $f;
}


sub drawBoxDown {
  my ($self, $f) = @_;
  $self->{DRAWBOX_DO_T} = $f;
}

#sub drawHorizontalScale {
#  my ($self, $f) = @_;
#  $self->{DRAWHORIZSCALE} = $f;
#}


sub setAlgoClust {
  my ($self, $f) = @_;
  $self->{ALGOCLUST} = $f;
}


sub setNeg {

  my ($self, $f) = @_;
  $self->{NEG} = $f;

}

sub drawFrames {
  my ($self, $f) = @_;
  $self->{DRAWFRAMES} = $f;
}

sub drawRowNames {
  my ($self, $f) = @_;
  $self->{DRAWROWNAMES} = $f;
}

sub colDrawMotifs {
  my ($self, $f) = @_;
  $self->{COLDRAWMOTIFS} = $f;
}

sub rowDrawMotifs {
  my ($self, $f) = @_;
  $self->{ROWDRAWMOTIFS} = $f;
}


sub addColumnDesc {
  my ($self, $file) = @_;

  my $ta = Table->new;
  $ta->loadFile($file);
  $self->{COLDESC} = $ta->getIndexKV(0,1);
}



sub addRowDesc {
  my ($self, $file) = @_;

  my $ta = Table->new;
  $ta->loadFile($file);
  $self->{ROWDESC} = $ta->getIndexKV(0,1);
}



sub setHeaderFontSize {
  my ($self, $f) = @_;
  $self->{HFONTSIZE}  = $f;
}

sub setRownameFontSize {
  my ($self, $f) = @_;
  $self->{RFONTSIZE}  = $f;
}


sub setXbase {
  my ($self, $f) = @_;
  $self->{XBASE}  = $f;
}

sub setXright {
  my ($self, $f) = @_;
  $self->{XRIGHT}  = $f;
}


sub setYbase {
  my ($self, $f) = @_;
  $self->{YBASE}  = $f;
}


sub setOutputFileName {
  my ($self, $f) = @_;
  $self->{OUTFILE}  = $f;
}

sub setH {
  my ($self, $m) = @_;
  $self->{H} = $m;
}

sub setW {
  my ($self, $m) = @_;
  $self->{W} = $m;
}


sub setVerbose {
  my ($self, $m) = @_;
  $self->{VERBOSE} = $m;
}


sub setMin {
  my ($self, $m) = @_;
  $self->{MIN} = $m;
}

sub setMax {
  my ($self, $m) = @_;
  $self->{MAX} = $m;
}

sub getMatrixDimensions {
  my ($self) = @_;
  my $ng     = scalar( @{$self->{MATRIX}} );
  my $nc     = scalar( @{$self->{MATRIX}->[0] } );
  
  return [ $ng, $nc ];

}


sub setAutoMinMax {
  my ($self) = @_;

  if (!defined($self->{MATRIX})) {
    die "Matrix undefined .. \n";
  }

  $self->{MIN} = $self->{MATRIX}->[0]->[0];
  $self->{MAX} = $self->{MATRIX}->[0]->[0];
  
  my $ng     = scalar( @{$self->{MATRIX}} );
  my $nc     = scalar( @{$self->{MATRIX}->[0] } );

   for (my $i=0; $i<$ng; $i++) {

    for (my $j=0; $j<$nc; $j++) {

      my $c = $self->{MATRIX}->[$i]->[$j];
      next if ( $c eq "");

      if ($c > $self->{MAX}) {
	$self->{MAX} = $c;
      }

      if ($c < $self->{MIN}) {
	$self->{MIN} = $c;
      }


    }
  }

}



sub setMinText {
  my ($self, $m) = @_;
  $self->{MINTEXT} = $m;
}

sub setMaxText {
  my ($self, $m) = @_;
  $self->{MAXTEXT} = $m;
}



sub pdfify {
  my ($self) = @_;
  my $feps = $self->{OUTFILE};
  my $fpdf = $feps;
  $fpdf =~ s/\.eps/\.pdf/;

  unlink $fpdf if ( -e $fpdf );

  my $todo = "ps2pdf  -dEPSCrop -dAutoRotatePages=/None  $self->{OUTFILE} $fpdf";
  system($todo);

  if (! -e $fpdf) {
    print "Warning: no PDF was created.\n";
  }
  
}


sub loadMatrix {
  my ($self, $filename, $header, $rownames) = @_;


  #if (!defined()) {
  #  $self->setColMap;
  #}


  my $ta = Table->new;
  $ta->loadFile($filename);

  my $a_ref = $ta->getArray();
  
  # if no header and no row 
  if ( ( ($header == 0) || !defined($header) ) && ( ($rownames == 0) || !defined($rownames) ) )   {


    # header and no rows
  } elsif ( ($header == 1) && ( ($rownames == 0) || !defined($rownames) ) )   {
    $self->{HEADER} = shift @$a_ref;


    # no header but rows
  } elsif ( ( ($header == 0) || !defined($header) ) && ($rownames == 1) )   {

    #print "NUH\n";

    $self->{ROWNAMES} = [];
    foreach my $r (@$a_ref) {
      my $n = shift @$r;
      push @{$self->{ROWNAMES}}, $n;
    }

    # header and row
  } else {

    $self->{HEADER} = shift @$a_ref;
    shift @{ $self->{HEADER} };

    $self->{ROWNAMES} = [];
    foreach my $r (@$a_ref) {
      my $n = shift @$r;
      push @{$self->{ROWNAMES}}, $n;
    }
  }

  $self->{MATRIX} = $a_ref;

  if ($self->{NEG} == 1) {
    
    foreach my $r (@{$self->{MATRIX}}) {
      foreach my $s (@$r) {
	$s = -1.0 * $s;
      }
    }
   
  } # neg

  if (defined($self->{COLTOCHANGESIGN})) {
    
    my $colidx = -1;
    for (my $i=0; $i<@{ $self->{HEADER}}; $i++) {
      if ($self->{HEADER}->[$i] eq $self->{COLTOCHANGESIGN}) {
	$colidx = $i;
      }
    }

    print "$colidx\n";

    for (my $i=0; $i<@{$self->{MATRIX}}; $i++) {
      $self->{MATRIX}->[$i]->[$colidx] = -1 * $self->{MATRIX}->[$i]->[$colidx];
    }
    
  }

}



sub loadRightTable {
  my ($self, $filename) = @_;



  my $ta = Table->new;
  $ta->loadFile($filename);

  $self->{RIGHTTABLE} = $ta->getArray();

}


sub setMatrix {
  my ($self, $a_ref, $r, $h) = @_;

  $self->{MATRIX}   = $a_ref;
  $self->{ROWNAMES} = $r;
  $self->{HEADER}   = $h;

}

sub regexpColNames {
  my ($self, $regexp) = @_;

  for (my $i=0; $i<@{$self->{HEADER}}; $i++) {
    my $txt = "\$self->{HEADER}->[\$i] =~ $regexp;";
    #print "$txt\n";
    eval($txt);

    #print "$self->{HEADER}->[$i]\n";

  }
  
}



sub regexpRowNames {
  my ($self, $regexp) = @_;

  for (my $i=0; $i<@{$self->{ROWNAMES}}; $i++) {

    my $txt = "\$self->{ROWNAMES}->[\$i] =~ $regexp;";
    eval($txt);

    if (defined($self->{ROWDESC})) {
      $txt = "\$self->{ROWDESC}->{ \$self->{ROWNAMES}->[\$i] } =~ $regexp;";
      eval($txt);
    }

  }
  
}


sub loadIclustPartition {
  my ($self, $filename) = @_;
  my $ta = Table->new;
  $ta->loadFile($filename);
  $ta->sortbycol(1);
  $self->{PARTITION} = $ta->getArray();
}


#
# use exprssion of 1 gene 
#
sub sortColumnsByGene {
  my ($self, $gene) = @_;

  my $a_ref_tmp = [];  # temporary matrix

  if ($self->{VERBOSE} == 1) {
    print "Reorganizing matrix cols using gene .. \n";
  }

  my $ng     = scalar( @{$self->{MATRIX}} );
  my $nc     = scalar( @{$self->{MATRIX}->[0] } );

  # get row
  my $idx = undef;
  for (my $i=0; $i<$ng; $i++) {
    print "$self->{ROWNAMES}->[$i]\n";
    if ($self->{ROWNAMES}->[$i] eq $gene) {
      $idx = $i;
      last;
    }
  }

  if (!defined($idx)) {
    print "$gene not found\n";
    return;
  }

  # make the order
  my $order = Sets::order($self->{MATRIX}->[$idx], 1);

  # sort matrix
  for (my $i=0; $i<$ng; $i++) {
    for (my $j=0; $j<$nc; $j++) {

      my $newj = $order->[$j];
      #print "newj=$newj\n";
      $a_ref_tmp->[$i]->[$j] = $self->{MATRIX}->[$i]->[$newj];

    }
  }

  $self->{MATRIX} = $a_ref_tmp;

  # sort header
  my $tmp = [];
  for (my $j=0; $j<$nc; $j++) {
    my $newj   = $order->[$j];
    $tmp->[$j] = $self->{HEADER}->[$newj];
  }
  $self->{HEADER} = $tmp;
}




#
# uses coldesc to reorder columns
#
sub reorderColumnsUsingDesc {
  my ($self) = @_;

  my $a_ref_tmp = [];  # temporary matrix

  if ($self->{VERBOSE} == 1) {
    print "Reorganizing matrix cols using coldesc.. \n";
  }

  my $ng     = scalar( @{$self->{MATRIX}} );
  my $nc     = scalar( @{$self->{MATRIX}->[0] } );

  # make the order
  #  1. get all labels, assign a unique id to each label
  my %H      = ();
  my $cnt    = 0;
  my @labels = values( %{ $self->{COLDESC} } );
  foreach my $l (@labels) {
    if (!defined($H{$l})) {
      $H{ $l } = $cnt;
      print "$l\n";
      $cnt++;
    }
  }	
  
  my @order = ();
  foreach my $lab (sort(keys(%H))) {
    for (my $i=0; $i<@{$self->{HEADER}}; $i++) {
      if ( $self->{COLDESC}->{ $self->{HEADER}->[$i] } eq $lab) {
	push @order, $i;
      }
    }
  }
  
      

  for (my $i=0; $i<$ng; $i++) {
    
    for (my $j=0; $j<$nc; $j++) {
      my $newj = $order[$j];
      
      $a_ref_tmp->[$i]->[$j] = $self->{MATRIX}->[$i]->[$newj];
    }
    
  }  
  
  $self->{MATRIX} = $a_ref_tmp;
  
  my $tmp = [];
  for (my $j=0; $j<$nc; $j++) {
    my $newj   = $order[$j];
    $tmp->[$j] = $self->{HEADER}->[$newj];
  }
  $self->{HEADER} = $tmp;
}


sub reorderColsUsingColList {

  my ($self, $colfile) = @_;

  my $a_ref_cols = Sets::readSet($colfile);
  
  my @order = ();
  my @newh  = ();
  foreach my $lab (@$a_ref_cols) {
    for (my $i=0; $i<@{$self->{HEADER}}; $i++) {
      if ( $self->{HEADER}->[$i] eq $lab) {
	push @order, $i;
	push @newh, $self->{HEADER}->[$i];
      }
    }
  }
  $self->{HEADER} = \@newh;

  print join("-", @order) . "\n";

  

  # update matrix  
  my $a_ref_tmp = [];  # temporary matrix

  my $ng     = scalar( @{$self->{MATRIX}} );
  my $nc     = scalar( @{$self->{MATRIX}->[0] } );

  for (my $i=0; $i<$ng; $i++) {
    foreach my $o (@order) {
      push @{ $a_ref_tmp->[$i] }, $self->{MATRIX}->[$i]->[$o];
    }
  }

  $self->{MATRIX} = $a_ref_tmp;

}


sub reorderUsingPartition {
  my ($self, $h, $v) = @_;

  my $a_ref_tmp = [];  # temporary matrix

  if ($self->{VERBOSE} == 1) {
    print "Reorganizing matrix .. \n";
  }

  my $ng     = scalar( @{$self->{MATRIX}} );
  my $nc     = scalar( @{$self->{MATRIX}->[0] } );
 
  for (my $i=0; $i<$ng; $i++) {
    my $newi = $self->{PARTITION}->[$i]->[0];

    for (my $j=0; $j<$nc; $j++) {
      my $newj = $self->{PARTITION}->[$j]->[0];
      
      $a_ref_tmp->[$i]->[$j] = $self->{MATRIX}->[$newi]->[$newj];
    }

  }  

  $self->{MATRIX} = $a_ref_tmp;
 
}


sub loadColumnsPartition {
  my ($self, $filename) = @_;

  my $ta = Table->new;
  $ta->loadFile($filename);

  # only keep the cluster numbers
  $self->{COLUMNS_PARTITION} = $ta->getColumn(1);
  shift @{ $self->{COLUMNS_PARTITION} };


}

sub clusterColumns {

  my ($self, $d) = @_;

  # transpose matrix to cluster rows
  my $a_ref_t = Sets::transpose($self->{MATRIX});


  my $a_ref_dist = [];
  my $ac = AggloClust->new;
  $ac->setAlgoClust($self->{ALGOCLUST});

  if ($d eq "euclidean") {
    $a_ref_dist = Sets::getPairwiseEuclideanDistanceMatrix($a_ref_t);
  } else {
    $a_ref_dist = Sets::getPairwisePearsonDistanceMatrix($a_ref_t);
  }

  $ac->setUseCorr(0);  
  $ac->setDistanceMatrix($a_ref_dist);
  
  $ac->agglomerate_using_max_linkage();

  my $a_ref_o = $ac->getDFSOrder();
  print join("-", @$a_ref_o) . "\n";

  my $n          = @{ $self->{MATRIX} };
  my $a_ref_newm = [];
  my $a_ref_newc = [];

  for (my $i=0; $i<$n; $i++) {
    for (my $j=0; $j<@$a_ref_o; $j++) {
      $a_ref_newm->[ $i ]->[ $j ] = $self->{MATRIX}->[ $i ]->[ $a_ref_o->[ $j ] ];
    }
  }

  for (my $j=0; $j<@$a_ref_o; $j++) {
    $a_ref_newc->[ $j ] = $self->{HEADER}->[ $a_ref_o->[ $j ] ];
  }

  $self->{MATRIX} = $a_ref_newm;
  $self->{HEADER} = $a_ref_newc;

}


sub clusterRows {

  my ($self, $d) = @_;

  my $a_ref_dist = undef;
  if ($d eq "euclidean") {
    $a_ref_dist = Sets::getPairwiseEuclideanDistanceMatrix($self->{MATRIX});
  } else {
    $a_ref_dist = Sets::getPairwisePearsonDistanceMatrix($self->{MATRIX});
  }

  my $ac = AggloClust->new;


  $ac->setDistanceMatrix($a_ref_dist);
  $ac->setUseCorr(0);
  #$ac->agglomerate_using_max_linkage();

  $ac->setAlgoClust($self->{ALGOCLUST});
  $ac->agglomerate_using_max_linkage();

  my $a_ref_o = $ac->getDFSOrder();
  my $n       = @{ $self->{MATRIX} };
  
  my $a_ref_newm = [];
  my $a_ref_newr = [];

  for (my $i=0; $i<@$a_ref_o; $i++) {
    $a_ref_newm->[ $i ] = $self->{MATRIX}->[ $a_ref_o->[ $i ] ];
    $a_ref_newr->[ $i ] = $self->{ROWNAMES}->[ $a_ref_o->[ $i ] ];
  }

  $self->{MATRIX}   = $a_ref_newm;
  $self->{ROWNAMES} = $a_ref_newr;

  # print join("-", @$a_ref_o) . "\n";

}



#
#
#
sub sortRowsByMin {
  my ($self, $rev) = @_;

  $rev = 0 if (!defined($rev));
    
  
  # number of rows
  my $n        = @{ $self->{MATRIX} };

    
  my @pe = ();
  for (my $i=0; $i<$n; $i++) {
    my $idx = Sets::indexMinInArray($self->{MATRIX}->[$i]);
    push @pe, $idx;
    #print "$idx\n";
  }

  my $o = Sets::order(\@pe, $rev);
    
  my $a_ref_newm = [];
  my $a_ref_newr = [];
    
  for (my $i=0; $i<@$o; $i++) {
    $a_ref_newm->[ $i ] = $self->{MATRIX}->  [ $o->[ $i ] ];
    $a_ref_newr->[ $i ] = $self->{ROWNAMES}->[ $o->[ $i ] ];
  }
  
  $self->{MATRIX}   = $a_ref_newm;
  $self->{ROWNAMES} = $a_ref_newr;
  
}



sub sortRowsByMax {
  my ($self, $rev) = @_;

  $rev = 0 if (!defined($rev));
    
  
  # number of rows
  my $n        = @{ $self->{MATRIX} };

    
  my @pe = ();
  for (my $i=0; $i<$n; $i++) {
    my $idx = Sets::indexMaxInSmoothedArray($self->{MATRIX}->[$i]);
    push @pe, $idx;
    print "$idx\n";
  }

  my $o = Sets::order(\@pe, $rev);
    
  my $a_ref_newm = [];
  my $a_ref_newr = [];
    
  for (my $i=0; $i<@$o; $i++) {
    $a_ref_newm->[ $i ] = $self->{MATRIX}->  [ $o->[ $i ] ];
    $a_ref_newr->[ $i ] = $self->{ROWNAMES}->[ $o->[ $i ] ];
  }
  
  $self->{MATRIX}   = $a_ref_newm;
  $self->{ROWNAMES} = $a_ref_newr;
  
}

#
# sort rows according to file
#
sub sortRowsUsingFile {
  my ($self, $file) = @_;

  # read file
  my $a_ref = [];
  open IN, $file or die "Cannot open $file\n";
  while (my $l = <IN>) {
    chomp $l;
    my @a = split /\t/, $l, -1;
    push @$a_ref, $a[0] if ($a[0] ne "");
  }
  close IN;

  
  # number of rows
  my $n        = @{ $self->{MATRIX} };

  my %POS = ();
  my $idx = 0;
  foreach my $r (@{$self->{ROWNAMES}}) {
    $POS{$r} = $idx;
    $idx ++;
  }

  # create order
  my $o = [];
  foreach my $m (@$a_ref) {
    push @$o, $POS{$m} if (defined($POS{$m}));
  }

  # make new matrix
  my $a_ref_newm = [];
  my $a_ref_newr = [];

  for (my $i=0; $i<@$o; $i++) {
    $a_ref_newm->[ $i ] = $self->{MATRIX}->  [ $o->[ $i ] ];
    $a_ref_newr->[ $i ] = $self->{ROWNAMES}->[ $o->[ $i ] ];
  }

  $self->{MATRIX}   = $a_ref_newm;
  $self->{ROWNAMES} = $a_ref_newr;

}



sub sortRowsUsingCorrelationWithGene {
  my ($self, $gene) = @_;

  # find gene
  my $n        = @{ $self->{MATRIX} };
  my $idx_gene = -1;

  for (my $i=0; $i<$n; $i++) {
    if ($self->{ROWNAMES}->[$i] eq $gene) {
      $idx_gene = $i;
      last;
    }
  }

  if ($idx_gene != -1) {

    
    my @pe = ();
    for (my $i=0; $i<$n; $i++) {
      my $c = Sets::pearson($self->{MATRIX}->[$i], $self->{MATRIX}->[$idx_gene]);

      push @pe, $c;
    }

    my $o = Sets::order(\@pe);
    
    my $a_ref_newm = [];
    my $a_ref_newr = [];
    
    for (my $i=0; $i<@$o; $i++) {
      $a_ref_newm->[ $i ] = $self->{MATRIX}->  [ $o->[ $i ] ];
      $a_ref_newr->[ $i ] = $self->{ROWNAMES}->[ $o->[ $i ] ];
    }
    
    $self->{MATRIX}   = $a_ref_newm;
    $self->{ROWNAMES} = $a_ref_newr;

    
  } else {
    die "Could not find $gene.\n";
  }
  

}


sub saveOrderedRowNames {
  my ($self, $f) = @_;

  Sets::writeSet($self->{ROWNAMES}, $f);

}

sub saveOrderedColNames {
  my ($self, $f) = @_;

  Sets::writeSet($self->{HEADER}, $f);

}

sub saveMatrix {
  my ($self, $f) = @_;

  my $n = @{$self->{ROWNAMES}};
  my $m = @{$self->{HEADER}};

  open OUT, ">$f" or die "Cannot open $f.\n";

  print OUT "ID\t" . join("\t", @{$self->{HEADER}}) . "\n";
  
  my $i = 0;
  foreach my $r (@{$self->{MATRIX}}) {
    print OUT $self->{ROWNAMES}->[$i] . "\t" . join("\t", @$r) . "\n";;
    $i++;
  }

  
  close OUT;


}



sub setColMap {

  my ($self, $cmapfile) = @_;

  my $ta = Table->new;
  $ta->setDelim("[\t\ ]");
  $ta->loadFile($cmapfile);
  $self->{COLMAP} = $ta->getArray();

  #if ($self->{VERBOSE} == 1) {
  #  foreach my $r (@{ $self->{COLMAP} }) {
  #    print "$r->[0]\n";
  #  }
  #}
  
  $ta->setDelim('\t');

}

sub drawHeader {
  my ($self, $n) = @_;
  $self->{DRAWHEADER} = $n;
}

sub drawHeaderArrows {
  my ($self, $txtl, $txtc, $txtr) = @_;
  
  $self->{ARROWS} = [ $txtl, $txtc, $txtr ];
  
}

sub draw {
  my ($self) = @_;

  if (@{ $self->{COLMAP} } == 0) {
    die "Colmap empty ... cannot continue, dying.\n";
  }

  my $ng     = scalar( @{$self->{MATRIX}} );
  my $nc     = scalar( @{$self->{MATRIX}->[0] } );
 
  if ($self->{VERBOSE} == 1) {
    print "Drawing a $ng x $nc matrix.\n";
  }	

  $self->{XSIZE}  = $self->{XBASE} + $self->{W} * $nc + $self->{XRIGHT};
  $self->{YSIZE}  = $self->{YBASE} + $self->{H} * $ng + 250;

  

  $self->{P} = new PostScript::Simple(xsize     => $self->{XSIZE},
				      ysize     => $self->{YSIZE},
 				      colour    => 1,
				      eps       => 1,
				      units     => "pt");
  
  $self->{P}->setlinewidth(0.25);

  my $maxlenrownames = 0;
  
  if (defined($self->{ROWDESC})) {
    # detrmine max rowname
    for (my $i=0; $i<@{$self->{ROWNAMES}}; $i++) {
      if (length($self->{ROWNAMES}->[$i]) > $maxlenrownames) {
	$maxlenrownames = length($self->{ROWNAMES}->[$i]);
      }
    }
    
  }

  #
  # define thresholds for density or log10 p-values
  #

  my @last_col = (-1, -1, -1);



  for (my $i=0; $i<$ng; $i++) {


    for (my $j=0; $j<$nc; $j++) {

      my $c = $self->{MATRIX}->[$i]->[$j];

      my @col = ();
      

      #print "$c\t$self->{ONLYSIGNIF}\t$self->{DRAWBOX_UP_T}'\n";

      if (($self->{ONLYSIGNIF} == 1) && (abs($c) < abs($self->{DRAWBOX_UP_T}))) {
	@col = "black";
	
      } elsif (($self->{ONLYSIGNIF} == 2) && ( (abs($c) < abs($self->{DRAWBOX_UP_T})) || ($c < 0)) ) {
	@col = "black";
	
      } elsif (($c eq "") || ($c eq "NA")) {
	@col = (127, 127, 127);

      } else {	
	if (!defined($self->{COLMAP})) {
	  @col = Sets::interp_general( $c, [0, 0, 204], [255, 255, 0], $self->{MIN}, $self->{MAX});
	} else {
	  @col = Sets::interp_from_matlab_colormap( $c, $self->{COLMAP}, $self->{MIN}, $self->{MAX});
	}

      } 


      #print "$self->{HEADER}->[$j]\tc=$c\t@col\n";


      if ( ($col[0] != $last_col[0]) ||
	   ($col[1] != $last_col[1]) ||
	   ($col[2] != $last_col[2]) ) {
	$self->{P}->setcolour(@col);
	@last_col = @col;
      }

      $self->{P}->box({ filled => 1},
		       $self->{XBASE} +  $j*$self->{W},
		       $self->{YSIZE} - ($self->{YBASE} +  $i*$self->{H}),
		       $self->{XBASE} +  $j*$self->{W}+$self->{W},
		       $self->{YSIZE} - ($self->{YBASE} +  $i*$self->{H} + $self->{H})
	      );
      


      #$self->{P}->box({filled => 1}, 
      #$self->{XBASE} + $j * $self->{W},      $self->{YSIZE} - ($self->{YBASE} + $i*$self->{H}) , 
      #	    $self->{XBASE} + $j * $self->{W} + $self->{W}, $self->{YSIZE} - ($self->{YBASE} + ($i*$self->{H}+$self->{H})));

    } # for (my $j
    
    if ($self->{VERBOSE} == 1) {
      my $perc = int( 0.5+100 * $i / $ng) + 1;
      if ($perc % 10 == 0) {
	print "$perc% drawn       \n";
      }
    }

    #last if ($i == 100);
    
  } # for (my $i


  
  if ($self->{DRAWSIGNIFBOXES} == 1) {

    $self->{P}->setlinewidth(0.75);
    for (my $i=0; $i<$ng; $i++) {
      for (my $j=0; $j<$nc; $j++) {
	
	my $c = $self->{MATRIX}->[$i]->[$j];
	
	my $draw = 0;

	if (($self->{DRAWSIGNIFBOXES} == 1) && ($c > 0) && (abs($c) > abs($self->{DRAWBOX_UP_T}))) {
	  $self->{P}->setcolour("red");
	  $draw = 1;
	}

	if (($self->{DRAWSIGNIFBOXES} == 1) && ($c < 0) && (abs($c) > abs($self->{DRAWBOX_UP_T}))) {
	  $self->{P}->setcolour("blue");
	  $draw = 1;
	}

	if ($draw == 1) {
	  $self->{P}->box({ filled => 0},
			  $self->{XBASE} +  $j*$self->{W},
			  $self->{YSIZE} - ($self->{YBASE} +  $i*$self->{H}),
			  $self->{XBASE} +  $j*$self->{W}+$self->{W},
			  $self->{YSIZE} - ($self->{YBASE} +  $i*$self->{H} + $self->{H})
			 );
	} # end if
      } # end for $j      
    } # end for $i
    
  } # end IF DRAWBOX


  if ($self->{DRAWFRAMES} == 1) {

    $self->{P}->setcolour("black");
    for (my $i=0; $i<$ng; $i++) {
      for (my $j=0; $j<$nc; $j++) {

	$self->{P}->box({ filled => 0},
			$self->{XBASE} +  $j*$self->{W},
			$self->{YSIZE} - ($self->{YBASE} +  $i*$self->{H}),
			$self->{XBASE} +  $j*$self->{W}+$self->{W},
			$self->{YSIZE} - ($self->{YBASE} +  $i*$self->{H} + $self->{H})
		       );
      }

    }

  } # for (my $i




  #
  # arrows
  #

  if (defined($self->{ARROWS})) {
    
    $self->{P}->setcolour("black");

    my $x1 = $self->{XBASE};
    my $x2 = $x1 + $nc * $self->{W};

    my $y1 = $self->{YSIZE} - ($self->{YBASE} - 10);
    my $y2 = $y1;

    $self->{P}->line($x1, $y1, $x2, $y2);

    # arrow left
    $self->{P}->line($x1, $y1, $x1+5, $y1-5);
    $self->{P}->line($x1, $y1, $x1+5, $y1+5);

    # arrow right
    $self->{P}->line($x2, $y2, $x2-5, $y2-5);
    $self->{P}->line($x2, $y2, $x2-5, $y2+5);


    $self->{P}->setfont($self->{FONT}, $self->{HFONTSIZE});
    $self->{P}->text({ align => "left"}, $x1+5, $y1+2, $self->{ARROWS}->[0]);
    $self->{P}->text({ align => "right" }, $x2-5, $y2+2, $self->{ARROWS}->[2]);

    $self->{P}->text({ align => "centre" }, ($x2+$x1)/2, $y2+20, $self->{ARROWS}->[1]);

  }

  #
  # DRAW HEADER
  #  
  if (defined($self->{HEADER}) && ($self->{DRAWHEADER} == 1)) {

    $self->{P}->setcolour("black");
    $self->{P}->setfont($self->{FONT}, $self->{HFONTSIZE});

    my $avg_header_len = 0;
    my $max_header_len = 0;

    for (my $j=0; $j<$nc; $j++) {
      $avg_header_len += length($self->{HEADER}->[$j]);
      if (length($self->{HEADER}->[$j]) > $max_header_len) {
	$max_header_len = length($self->{HEADER}->[$j]);
      }
    }
    $avg_header_len /= @{$self->{HEADER}};

    for (my $j=0; $j<$nc; $j++) {


      if ($self->{COLDRAWMOTIFS} == 0) {

	my $h = $self->{HEADER}->[$j];

	if (defined($self->{COLDESC})) {
	  
	  my $hh = $self->{COLDESC}->{ $self->{HEADER}->[$j] };
	  
	  $hh =~ s/\.txt//;
	  $hh =~ s/M\d{5}\_//;
	  #$self->{P}->text({ rotate => $self->{COLNAME_ANGLE} },
	  #$self->{XBASE} +  $j * $self->{W} + 3 * $self->{W} / 4,
	  #		   $self->{YSIZE} -  ( $self->{YBASE} - 3 - $max_header_len * (0.7) * $self->{HFONTSIZE}) ,
	  #		   $hh 
	  #);
	  $h = "$hh";
	}	


	$self->{P}->text({ rotate => $self->{COLNAME_ANGLE} },
			 $self->{XBASE} +  $j * $self->{W} + $self->{W} / 2 + $self->{HFONTSIZE}*(1/4),
			 $self->{YSIZE} -  ( $self->{YBASE} - 3),
			 $h
			);
	

	
	
      } else {
	
	if (! -e "TMP") {
	  mkdir "TMP";
	}
	
	my $mo = Sets::myre2wm($self->{HEADER}->[$j]);
	
	my $prefix = "TMP/";
	my $tmptxt = "$prefix$j.txt";
	my $tmpeps = "$prefix$j.eps";

	open OUT, ">$tmptxt" or die "cannot open $tmptxt\n";
	print OUT $mo;
	close OUT;

	if ($self->{COLDRAWMOTIFS} == 1) {
	  system("$self->{WEBLOGODIR}/weblogo/seqlogo -f $tmptxt -F EPS  -a -c -M -n -Y -w 5 -h 3 > $tmpeps");
	}
	
	my $e  = new PostScript::Simple::EPS(file => $tmpeps);

	# get height
	my $eh = $e->height;
	my $ew = $e->width;
	
	my $fa = 1.25;

	# height must be $h, so scale down to $h = k * LO
	$e->scale( ($self->{W}*$fa) / $eh);
	$e->rotate(90);
	my $ew_new = int(0.5 + $ew * $self->{W} * $fa / $eh);

	# finally add to picture
	$self->{P}->_add_eps($e, $self->{XBASE} + $self->{W} * $fa + $j * $self->{W},  $self->{YSIZE} - $self->{YBASE} ); 

	

	if (defined($self->{COLDESC})) {

	  my $hh = $self->{COLDESC}->{ $self->{HEADER}->[$j] };

	  $hh =~ s/\.txt//;
	  $hh =~ s/M\d{5}\_//;
	  $self->{P}->text({ rotate => 90 },
			   $self->{XBASE} +  $j * $self->{W} + 3 * $self->{W} / 4,
			   $self->{YSIZE} -  ( $self->{YBASE} - 3 - $ew_new ),
			   $hh 
			  );

	}	

	
      }
	
    }
    
  }
  
  
  #
  # DRAW ROWNAMES
  #

  if (defined($self->{ROWNAMES})) {
    
    $self->{P}->setcolour("black");
    $self->{P}->setfont($self->{FONT}, $self->{RFONTSIZE});  # 5
    
    for (my $i=0; $i<$ng; $i++) {
      
      #print "Row name $self->{ROWNAMES}->[$i]\n";

      #$self->{ROWNAMES}->[$i] =~ s/\ .+$//;

      

      if ($self->{DRAWROWNAMES} eq 'left') {
	$self->{P}->text({align => 'right'},,
			 $self->{XBASE} - 3,  
			 $self->{YSIZE} -  ( $self->{YBASE} + $i * $self->{H} + int(0.5+$self->{H}/2)),
			 $self->{ROWNAMES}->[$i]
			);
      } else {

	my $ew_new = undef;
	if ($self->{ROWDRAWMOTIFS} == 0) {
	
	  $self->{P}->text({align => 'left'},,
			   $self->{XBASE} + $nc * $self->{W} + 3,  
			   $self->{YSIZE} -  ( $self->{YBASE} + $i * $self->{H} + int(0.5+$self->{H}/2) +  $self->{RFONTSIZE}/2.5),
			   $self->{ROWNAMES}->[$i]
			  );
	  
	} else {
	  
	  # draw motifs in columns
	  
	  if (! -e "TMP") {
	    mkdir "TMP";
	  }
	
	  my $mo = Sets::myre2wm($self->{ROWNAMES}->[$i]);
	
	  my $prefix = "TMP/";
	  my $tmptxt = "$prefix$i.txt";
	  my $tmpeps = "$prefix$i.eps";
	  
	  open OUT, ">$tmptxt" or die "cannot open $tmptxt\n";
	  print OUT $mo;
	  close OUT;
	  
	  if ($self->{ROWDRAWMOTIFS} == 1) {
	    system("$self->{WEBLOGODIR}/weblogo/seqlogo -f $tmptxt -F EPS  -a -c -M -n -Y -w 5 -h 3 > $tmpeps");
	  }
	
	  my $e  = new PostScript::Simple::EPS(file => $tmpeps);
	  
	  # get height
	  my $eh = $e->height;
	  my $ew = $e->width;
	  
	  my $fa = 1.25;
	  
	  # height must be $h, so scale down to $h = k * LO
	  $e->scale( ($self->{H}*$fa) / $eh);
	  #$e->rotate(90);
	  $ew_new = int(0.5 + $ew * $self->{H} * $fa / $eh);
	  
	  # finally add to picture
	  $self->{P}->_add_eps($e, $self->{XBASE} + $nc * $self->{W},  $self->{YSIZE} - ( $self->{YBASE} + $i * $self->{H}  + $self->{H} + 0.25 * $self->{H}) ); 

	}

	
	if (defined($self->{ROWDESC})) {

	  my $rd = $self->{ROWDESC}->{ $self->{ROWNAMES}->[$i] };

	  #my $nsp = $maxlenrownames - length($self->{ROWNAMES}->[$i]);	  
	  #print "MAX = $maxlenrownames, NSP=$nsp\n";
	  #my $spaces = " " x $nsp;
	  #$rd = $spaces . $rd;	  
	  #print "\"$rd\"\n";

	  my $x  = undef; 
	  if ($self->{ROWDRAWMOTIFS} == 0) {
	    $x = $self->{XBASE} + $nc * $self->{W} + 10 + $maxlenrownames *(3/4)*$self->{HFONTSIZE};
	  } else {
	    $x = $self->{XBASE} + $nc * $self->{W} + 10 + $ew_new; 
	  }

	  
	  my $y = $self->{YSIZE} -  ( $self->{YBASE} + $i * $self->{H} + int(0.5+$self->{H}/2) + 6);
	  

	  
	  $self->{P}->text({align => 'left'},
			   $x,
			   $y,
			   $rd
			);
	
	  
	}
	
      }
      
      
      
    }

  } # draw rownames


  
  if (defined($self->{RIGHTTABLE})) {
    
    $self->{P}->setfont($self->{FONT}, $self->{RFONTSIZE});  # 5

    my $tabh  = shift @{$self->{RIGHTTABLE}};
    
    my $ncols = @{$self->{RIGHTTABLE}->[0]};
    my $nrows = @{$self->{RIGHTTABLE}};

    my $offsetstart = 3 + 30;
    my $colwidth    = 40;

    for (my $j=0; $j<$ncols; $j++) {
      $self->{P}->setcolour("black");
      $self->{P}->text({align => 'left'},,
		       $self->{XBASE} + $nc * $self->{W} + $offsetstart + $j * $colwidth,  
		       $self->{YSIZE} - ( $self->{YBASE} - 3),
		       $tabh->[$j]
		      );
    }
    

    for (my $i=0; $i<$nrows; $i++) {
      
      if ($self->{RIGHTTABLE}->[$i]->[0] > 0) {
	$self->{P}->setcolour("red");	
      } else {
	$self->{P}->setcolour("green");
      }
      
      $self->{P}->box({ filled => 1},
		      $self->{XBASE} +  $nc*$self->{W} + $offsetstart - 5,
		      $self->{YSIZE} - ($self->{YBASE} +  $i*$self->{H}),
		      $self->{XBASE} +   $nc*$self->{W} + $offsetstart + $colwidth * $ncols,
		      $self->{YSIZE} - ($self->{YBASE} +  $i*$self->{H} + $self->{H})
		     );

      for (my $j=0; $j<$ncols; $j++) {

	$self->{P}->setcolour("black");
	$self->{P}->text({align => 'left'},,
			 $self->{XBASE} + $nc * $self->{W} + $offsetstart + $j * $colwidth,  
			 $self->{YSIZE} -  ( $self->{YBASE} + $i * $self->{H} + int(0.5+$self->{H}/2) +  $self->{RFONTSIZE}/2.5),
			 $self->{RIGHTTABLE}->[$i]->[$j]
			);
      }
      
    }
  } # right table
       
}


sub addClusters {
  my ($self, $ho, $ve) = @_;

  my $ng     = scalar( @{$self->{MATRIX}} );
  my $nc     = scalar( @{$self->{MATRIX}->[0] } );

  for (my $i=1; $i<$ng; $i++) {

    if ($self->{CLUSTER}->[$i] != $self->{CLUSTER}->[$i-1]) {
      
      $self->{P}->setcolour('white');
      $self->{P}->line(
		       $self->{XBASE},
		       $self->{YSIZE} - ($self->{YBASE} + $i  * $self->{H}),
		       $self->{XBASE} + $nc * $self->{W},
		       $self->{YSIZE} - ($self->{YBASE} + $i  * $self->{H})
		      );
      
      #
      # add Cx
      #
      # my $s = $CSIZES{$pref_c};    
      # calculate where it should be:
      # my $y = $i - int( 0.5 + $s / 2 ) + 5;
      # $self->{P}->setcolour("black");
      # $self->{P}->setfont("Courrier", 18);
      # $self->{P}->text({align => 'right'}, $self->{XBASE}, $self->{YSIZE} - ($self->{YBASE} + $y), "C$pref_c");
      #
      
    }
    
    for (my $i=1; $i<$nc; $i++) {
      
      if ($self->{CLUSTER}->[$i] != $self->{CLUSTER}->[$i-1]) {
	$self->{P}->line(
			 $self->{XBASE} + $i  * $self->{W}, 
			 $self->{YSIZE} - ($self->{YBASE} + 0),
			 $self->{XBASE} + $i  * $self->{W}, 
		       $self->{YSIZE} - ($self->{YBASE} + $ng * $self->{H})
			); 
      }
    }
  }
}


sub addColumnClusters {
  my ($self) = @_;

  my $ng     = scalar( @{$self->{MATRIX}} );
  my $nc     = scalar( @{$self->{MATRIX}->[0] } );

  $self->{P}->setcolour('white');
  $self->{P}->setlinewidth(0.25);

  for (my $i=1; $i<$nc; $i++) {

    if ($self->{COLUMNS_PARTITION}->[$i] != $self->{COLUMNS_PARTITION}->[$i-1]) {
      $self->{P}->line(
	       $self->{XBASE} +  $i * $self->{W},
	       $self->{YSIZE} - ($self->{YBASE} + 0),
	       $self->{XBASE} +  $i * $self->{W},
	       $self->{YSIZE} - ($self->{YBASE} + $ng * $self->{H})
	    );
    }

  }



}



sub output {
  my ($self) = @_;
  $self->{P}->output($self->{OUTFILE});
}


sub drawScale {
  my ($self, $x, $y) = @_;
  
  my $h   = 2; 
  my $w   = 20;
  my $sep = 0;
  my $res = 50;

  $self->{P}->setcolour("black");
  $self->{P}->setfont($self->{FONT}, 14);
  $self->{P}->text({align => "center"}, $x+$w/2+0, $self->{YSIZE} - ($y - 3), $self->{MAX});
  $self->{P}->setfont($self->{FONT}, 14);
  $self->{P}->text({align => "left", rotate => 90}, $x+$w/2+5, $self->{YSIZE} - ($y - 17), $self->{MAXTEXT});


  my $t = $self->{MAX};

  for (my $i=0; $i<=$res; $i++) {
    
    my @col = ();
    if (!defined($self->{COLMAP})) {
      @col = Sets::interp_general( $t, [0, 0, 204], [255, 255, 0], $self->{MIN}, $self->{MAX});
    } else {
      @col = Sets::interp_from_matlab_colormap( $t, $self->{COLMAP}, $self->{MIN}, $self->{MAX});
    }
    
    $self->{P}->setcolour( @col );
    $self->{P}->box({filled => 1}, $x, $self->{YSIZE} - ($y + $sep + $i*$h) , $x+$w, $self->{YSIZE} - ($y + $sep + $i*$h + $h));
   
    $t -= ($self->{MAX} - $self->{MIN}) / $res;
    
  }

  $self->{P}->setcolour( 'black' );
  $self->{P}->setfont($self->{FONT}, 14);

  $self->{P}->text({align => "center"}, $x+$w/2+1, $self->{YSIZE} - ($y + $sep + $res*$h + 14), $self->{MIN});

  

  $self->{P}->setfont($self->{FONT}, 14);
  $self->{P}->text({align => "right", rotate => 90}, $x+$w/2+5, $self->{YSIZE} - ($y + $sep + $res*$h + 11 + 10), $self->{MINTEXT});

  
}


#  drawHorizontalScale($xbase + @{$a_ref_M->[0]} * $w / 2 - 100 - 100, $ybase + @$a_ref_M*$h+ $h*8/4 , $min, $max, 100, $p, $xsize, $ysize, $colmap, $A_REF_COLMAP, "Enrichment", "Depletion", $scalefont, 20, 2, 1, 2);

sub drawHorizontalScale {
 my ($self, $x, $y) = @_;

 my $sep       = 0;
 my $scalefont = $self->{SCALEFONTSIZE};
 my $res       = $self->{SCALERES};
 my $h   = 20; 
 my $w   = 2;


 if ($self->{DRAWP} == 0) {
   # print MAX
   $self->{P}->setcolour("black");
   $self->{P}->setfont($self->{FONT}, $scalefont);
   $self->{P}->text({align => "right"}, $x-3, $self->{YSIZE} - ($y+$h-3), $self->{MAX});
   
   # print TEXT that goes with MAX
   $self->{P}->setfont($self->{FONT}, $scalefont);
   $self->{P}->text({align => "right"}, $x - 3 - $scalefont, $self->{YSIZE} - ($y+$h-3), $self->{MAXTEXT});

 } else {

   
   if (($self->{MIN} < 0) && ($self->{ONLYPOS} == 0)) {
     # print MAX
     $self->{P}->setcolour("black");
     $self->{P}->setfont($self->{FONT}, $scalefont);
     #$self->{P}->text({align => "right"}, $x-5, $self->{YSIZE} - ($y+$h-3), $self->{MAXTEXT});
     $self->{P}->text({align => "center"}, $x, $self->{YSIZE} - ($y+$h+$scalefont), $self->{MAXTEXT});
   }
 }

 my $t = $self->{MAX};

 my $fpv = 0; 

 # tmp
 #$self->{MINABSLP} = undef;

 for (my $i=0; $i<=$res; $i++) {

   my $dec = ($self->{MAX} - (defined($self->{MINABSLP})?$self->{MINABSLP}:$self->{MIN})) / $res;

  

   my @col = () ;
   if (!defined( $self->{COLMAP} )) {
     if ($i>$res/2) {
       @col = Sets::interp_general( $t, [0, 0, 0], [255, 0, 0], $self->{MIN}, 0);
     } else {
       @col = Sets::interp_general( $t, [255, 0, 0], [255, 255, 0], 0, $self->{MAX});
     }
   } else {
     @col = Sets::interp_from_matlab_colormap( $t, $self->{COLMAP}, $self->{MIN}, $self->{MAX});
   }

   $self->{P}->setcolour( @col );
   $self->{P}->box({filled => 1}, $x + $i*$w, $self->{YSIZE} - ($y) , $x + ($i+1)*$w  , $self->{YSIZE} - ($y + $h));

   $t -= $dec;

 }

 
 # second pass, for drawing p-values

 if ($self->{DRAWP} == 1) { 
   
   $t = $self->{MAX};
   $fpv = 0;

   my $min = (defined($self->{MINABSLP})?$self->{MINABSLP}:$self->{MIN});

   my $dec = ($self->{MAX} - $min) / $res;

   $self->{P}->setcolour( 'black' );
   my $smallfont = $scalefont*(3/4);
   $self->{P}->setfont($self->{FONT}, $smallfont);

   
   # positive p-values
   for (my $i=$self->{MINABSLP}; $i<=$self->{MAX}; $i+=$self->{LPINCREMENTS}) {
     
     my $txt = undef;
     if ($i == 0) {
       $txt = "p=1";	 
     } elsif ($i == $self->{MAX}) {
       $txt = "p<=1e-$i";
     } else {
       $txt = "p=1e-$i";
     }
     
     # total lentgh = $res * $w
     
     # 0 => a MAX + b
     
     # ?   =>  i
     
     # res w => a MIN + b
     
     my $a = $res * $w / ( $min - $self->{MAX});
     my $b = - $a * $self->{MAX};
     

     # 

     $self->{P}->text({align => "left", rotate => 45}, 
		      $x+$a*$i+$b+$smallfont/4,
		      $self->{YSIZE} - ($y - 2),
		      $txt);

     
   }


   if ($self->{ONLYPOS} == 0) {
     #
     # NEGATIVE p-values
     #
     for (my $i=-$self->{LPINCREMENTS}; $i>=$self->{MIN}; $i-=$self->{LPINCREMENTS}) {
       
       my $txt = undef;
       my $ii  = abs($i); 
       if ($i == 0) {
	 $txt = "p=1";	 
       } elsif ($i == $self->{MIN}) {
	 $txt = "p<=1e-$ii";
       } else {
	 $txt = "p=1e-$ii";
       }
       
       # total lentgh = $res * $w
       
       # 0 => a MAX + b
       
       # ?   =>  i
       
       # res w => a MIN + b
       
       my $a = $res * $w / ( $min - $self->{MAX});
       my $b = - $a * $self->{MAX};
       
       
       # 
       
       $self->{P}->text({align => "left", rotate => 45}, 
			$x+$a*$i+$b+$smallfont/4,
			$self->{YSIZE} - ($y - 2),
			$txt);
     }
   
     
   }

 }


 if ($self->{DRAWP} == 0) {
   # print MIN
   $self->{P}->setcolour("black");
   $self->{P}->setfont($self->{FONT}, $scalefont);
   $self->{P}->text({align => "left"}, $x+3+$res*$w,  $self->{YSIZE} - ($y+$h-3), $self->{MIN});
   
   # print TEXT that goes with MIN
   $self->{P}->setfont($self->{FONT}, $scalefont);
   $self->{P}->text({align => "left"}, $x + 3+$res*$w + $scalefont,  $self->{YSIZE} - ($y+$h-3), $self->{MINTEXT});

 } else {

   if (($self->{MIN} >= 0) || ($self->{ONLYPOS} == 1)) {
     $self->{P}->setcolour("black");
     $self->{P}->setfont($self->{FONT}, $scalefont);
     $self->{P}->text({align => "center"}, $x+$res*$w/2,  $self->{YSIZE} - ($y+$h+$scalefont), "Enrichment significance");
   } else {
     # print MIN
     $self->{P}->setcolour("black");
     $self->{P}->setfont($self->{FONT}, $scalefont);
     # WAS $self->{P}->text({align => "left"}, $x+5+$res*$w,  $self->{YSIZE} - ($y+$h-3), $self->{MINTEXT});
     $self->{P}->text({align => "center"}, $x+$res*$w,  $self->{YSIZE} - ($y+$h+$scalefont), $self->{MINTEXT});
   }
 


   # additional black box
   if (defined($self->{ONLYSIGNIF}) && ($self->{ONLYSIGNIF} == 1)) {
     
     $self->{P}->setcolour("black");
     my $xx = $x + $res*$w + $h*2;
     $self->{P}->box({filled => 1}, $xx,  $self->{YSIZE} - ($y) , $xx+$h,  $self->{YSIZE} - ($y + $h));

     $self->{P}->setcolour("black");
     $self->{P}->setfont($self->{FONT}, $scalefont);

     $self->{P}->text($xx + $h + 3, $self->{YSIZE} - ($y + $h - 3), "Non-significant"); 

   }

 }
   

}



1;

