package MyBarPlot;

use lib "$ENV{FIREDIR}/SCRIPTS/PostScript-Simple-0.07/lib";
use Table;
use Sets;
use PostScript::Simple;


sub new {
  my ($class) = @_;
  my $self    = {};

  $self->{MATRIX}   = undef;
  $self->{LABELS}   = undef;

  $self->{W}        = 20;
  $self->{MAXH}     = 50;

  $self->{INT_X}    = 20;
  $self->{INT_Y}    = 20;
  
  $self->{MIDLINES} = undef;

  $self->{COLMAP}   = undef;

  $self->{XBASE}    = 55;
  $self->{YBASE}    = 200;

  $self->{P}        = undef;
  $self->{XRIGHT}   = 50;

  $self->{VERBOSE}  = 1;
  $self->{OUTFILE}  = undef;
  
  $self->{MIN}      = 0;
  $self->{MAX}      = 1;


  $self->{COLOR}    = "black";

  $self->{COLDRAWMOTIFS}     = 1;
  $self->{WEBLOGODIR}        = "$ENV{HOME}/PROGRAMS/FIRE/SCRIPTS";

  $self->{ROWMOTIF} = 0;
  
  bless($self);
  return $self;
}

sub setRowMotif {
  my ($self, $m) = @_;
  $self->{ROWMOTIF} = $m;
}

sub setXRight {
  my ($self, $m) = @_;
  $self->{XRIGHT} = $m;
}

sub setMatrix {
  my ($self, $m) = @_;
  $self->{MATRIX} = $m;
}


sub setHeader {
  my ($self, $m) = @_;
  $self->{HEADER} = $m;
}

sub setRowNames {
  my ($self, $m) = @_;
  $self->{ROWNAMES} = $m;
}

sub setMinMax {
  my ($self, $min, $max) = @_;
  $self->{MAX} = $max;
  $self->{MIN} = $min;
}

sub setLabels {
  my ($self, $m) = @_;
  $self->{LABELS} = $m;
}


sub setBarColor {
  my ($self, $m) = @_;
  
  $self->{COLOR} = $m;
  
}

sub setOutputFileName {
  my ($self, $f) = @_;
  $self->{OUTFILE}  = $f;
}

sub setMidLines {
  my ($self, $m) = @_;
  $self->{MIDLINES} = $m;
}

#
# draw 
#
sub draw {
  my ($self) = @_;

  #if (@{ $self->{COLMAP} } == 0) {
  #  die "Colmap empty ... cannot continue, dying.\n";
  #}

  # number of rows and columns
  my $nr     = scalar( @{$self->{MATRIX}} );
  my $nc     = scalar( @{$self->{MATRIX}->[0] } );

  if ($self->{VERBOSE} == 1) {
    print "Drawing a $nr x $nc matrix.\n";
  }	

  $self->{XSIZE}  = $self->{XBASE} + ($self->{W}    + $self->{INT_X}) * $nc + $self->{XRIGHT};
  $self->{YSIZE}  = $self->{YBASE} + ($self->{MAXH} + $self->{INT_Y}) * $nr + 50;

  $self->{P} = new PostScript::Simple(xsize     => $self->{XSIZE},
				      ysize     => $self->{YSIZE},
 				      colour    => 1,
				      eps       => 1,
				      units     => "pt");

  $self->{P}->setlinewidth(1);


  #
  # header
  #
  for (my $j=0; $j<$nc; $j++) {
    
    my $x1 = $self->{XBASE} + $j * ($self->{W} + $self->{INT_X}) + $self->{W}/2 + 5;
    my $y1 = $self->{YSIZE} - ( $self->{YBASE} - 25 );
    
    $self->{P}->setcolour("black");
    $self->{P}->setfont("Times", 14);
    $self->{P}->text({ rotate => 90}, $x1, $y1, $self->{HEADER}->[$j] );
  } # for (my $j
  
  #
  # cycle thru rows
  #
  for (my $i=0; $i<$nr; $i++) {

    for (my $j=0; $j<$nc; $j++) {

      my $f  = $self->{MATRIX}->[$i]->[$j];

      my $h  = $self->{MAXH} * $f / $self->{MAX};

      my $x1 = $self->{XBASE} + $j * ($self->{W} + $self->{INT_X});
      my $x2 = $self->{XBASE} + $j * ($self->{W} + $self->{INT_X}) + $self->{W};
      
      my $y1 = $self->{YSIZE} - ( $self->{YBASE} + $i * ($self->{MAXH} + $self->{INT_Y}) + $self->{MAXH}  );
      my $y2 = $self->{YSIZE} - ( $self->{YBASE} + $i * ($self->{MAXH} + $self->{INT_Y}) + $self->{MAXH} - $h );

      $self->{P}->setcolour($self->{COLOR});
      $self->{P}->box({ filled => 1}, $x1, $y1, $x2, $y2 );

      if (defined($self->{LABELS})) {
	$self->{P}->setcolour("black");
	$self->{P}->setfont("Times", 8);
	my $txt = $self->{LABELS}->[$i]->[$j];
	$self->{P}->text({align => "centre"}, $x1 + $self->{W}/2, $y2+5, $txt);
      }

    } # for (my $j

    $self->{P}->setcolour("black");
   
    # bottom line
    my $x1 = $self->{XBASE} - $self->{W}/2; 
    my $x2 = $self->{XBASE} + $nc * ($self->{W} + $self->{INT_X}) - $self->{INT_X} + $self->{W}/2;      
    my $y1 = $self->{YSIZE} - ( $self->{YBASE} + $i * ($self->{MAXH} + $self->{INT_Y}) + $self->{MAXH}  );
    $self->{P}->line($x1, $y1, $x2, $y1);

    # intermediate line
    if (defined($self->{MIDLINES})) {
      
      my $h  = $self->{MAXH} * $self->{MIDLINES}->[$i] / $self->{MAX};	    
      my $y1 = $self->{YSIZE} - ( $self->{YBASE} + $i * ($self->{MAXH} + $self->{INT_Y}) + $self->{MAXH} - $h );
      $self->{P}->setcolour("black");
      $self->{P}->line($x1, $y1, $x2, $y1);

      $self->{P}->setcolour("black");
    }

    # side line
    my $x1 = $self->{XBASE} - $self->{W}/2; 
    my $y1 = $self->{YSIZE} - ( $self->{YBASE} + $i * ($self->{MAXH} + $self->{INT_Y}) + $self->{MAXH}  );
    my $y2 = $self->{YSIZE} - ( $self->{YBASE} + $i * ($self->{MAXH} + $self->{INT_Y}) + $self->{MAXH} - $self->{MAXH} );
    $self->{P}->line($x1, $y1, $x1, $y2);

    # upper tick line
    my $x1 = $self->{XBASE} - $self->{W}/2 - 3;
    my $x2 = $self->{XBASE} - $self->{W}/2 + 3;
    my $y1 = $self->{YSIZE} - ( $self->{YBASE} + $i * ($self->{MAXH} + $self->{INT_Y}) + $self->{MAXH}  - $self->{MAXH}  );
    $self->{P}->line($x1, $y1, $x2, $y1);
    $self->{P}->setcolour("black");
    $self->{P}->setfont("Courrier", 10);
    $self->{P}->text({ align => "right"}, $x1, $y1-3, $self->{MAX});

    # upper tick line
    my $x1 = $self->{XBASE} - $self->{W}/2 - 3;
    my $x2 = $self->{XBASE} - $self->{W}/2 + 3;
    my $y1 = $self->{YSIZE} - ( $self->{YBASE} + $i * ($self->{MAXH} + $self->{INT_Y}) + $self->{MAXH}   );
    $self->{P}->line($x1, $y1, $x2, $y1);
    $self->{P}->setcolour("black");
    $self->{P}->setfont("Courrier", 10);
    $self->{P}->text({ align => "right"}, $x1, $y1-3, $self->{MIN});

    # row names

    if (defined($self->{ROWNAMES})) {

      my $x2 = $self->{XBASE} + $nc * ($self->{W} + $self->{INT_X}) - $self->{INT_X} + $self->{W}/2;      
      my $y2 = $self->{YSIZE} - ( $self->{YBASE} + $i * ($self->{MAXH} + $self->{INT_Y}) + $self->{MAXH}  );

      if ($self->{ROWMOTIF} == 0) {	
	$self->{P}->setcolour("black");
	$self->{P}->setfont("Courrier", 14);
	$self->{P}->text($x2, $y2, $self->{ROWNAMES}->[$i]);

      } else {
	
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

	if ($self->{COLDRAWMOTIFS} == 1) {
	  system("$self->{WEBLOGODIR}/weblogo/seqlogo -f $tmptxt -F EPS  -a -c -M -n -Y -w 5 -h 3 > $tmpeps");
	}
	
	my $e  = new PostScript::Simple::EPS(file => $tmpeps);

	# get height
	my $eh = $e->height;
	my $ew = $e->width;
	
	my $fa = 1;

	# height must be $h, so scale down to $h = k * LO
	$e->scale( ($self->{MAXH}*$fa) / $eh);
	#$e->rotate(90);
	my $ew_new = int(0.5 + $ew * $self->{MAXH} * $fa / $eh);

	# finally add to picture
	$self->{P}->_add_eps($e, $x2, $y2 - ($self->{MAXH} - $ew_new) - $self->{MAXH} ); 

	

	
      }

    }

  } # for (my $i


}


sub output {
  my ($self) = @_;
  $self->{P}->output($self->{OUTFILE});
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



1;
