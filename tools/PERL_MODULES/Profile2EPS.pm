package Profile2EPS;


use lib "$ENV{FIREDIR}/SCRIPTS/PostScript-Simple-0.07/lib";
use Table;
use Sets;
use PostScript::Simple;
use strict;


sub new {
  
  my $self = {};

  $self->{PROFILES} = [];


  $self->{W}        = 200;
  $self->{H}        = 100;

  $self->{MARGIN}   = 10;
  
  $self->{XBASE}    = 50;
  $self->{YBASE}    = 50;
  $self->{XRIGHT}   = 50;
  $self->{ADD}      = 0;

  $self->{HLINE}    = undef;

  bless($self);
  return $self;
  
}


sub addProfile {
  my ($self, $r, $mult) = @_;
  
  my @a = @$r;

  if (defined($mult)) {
    foreach my $s (@a) {
      $s->[1] = $s->[1] * $mult;
    }
  }

  push @{ $self->{PROFILES} }, \@a;

  


}


sub addLine {
  my ($self, $h) = @_;

  $self->{HLINE} = $h;
  
}


sub interpolate {
  my ($self) = @_;
  
  my @b = ();
  foreach my $r (@{ $self->{PROFILES} }) {

    my $s = {};
    push @$s, $r->[0];
    for (my $i=1; $i<@$r; $i++) {
      my @a = ( 
	       ($r->[$i-1]->[0] + $r->[$i-1]->[0])/2,
	       ($r->[$i-1]->[1] + $r->[$i-1]->[1])/2
	      );
      push @$s, \@a;
      push @$s, $r->[$i];
    }	

    push @b, $s;
  }

  $self->{PROFILES} = \@b;
  
}


sub draw {

  my ($self) = @_;

  $self->{XSIZE}  = $self->{XBASE} + $self->{W} + $self->{XRIGHT};
  $self->{YSIZE}  = $self->{YBASE} + $self->{H} + 100;

  

  $self->{P} = new PostScript::Simple(xsize     => $self->{XSIZE},
				      ysize     => $self->{YSIZE},
 				      colour    => 1,
				      eps       => 1,
				      units     => "pt");



  
  my $n    = @{ $self->{PROFILES}->[0] };
  my $r    =    $self->{PROFILES}->[0];

  # determine min(X) and max(X)
  my $minX      = $r->[0]->[0];
  my $maxX      = $r->[$n-1]->[0];
  my $d_real    = $maxX - $minX;
  my $d_here    = $self->{W}; 
  my $coefX     = $d_here / $d_real;

  # same for Y
  my $minY = $r->[0]->[1];
  my $maxY = $r->[0]->[1];
  for (my $i=1; $i<$n; $i++) {
    if ($r->[$i]->[1] < $minY) {
      $minY = $r->[$i]->[1];
    }
    if ($r->[$i]->[1] > $maxY) {
      $maxY = $r->[$i]->[1];
    }
  }

  #print "$minY / $maxY\n";

  $maxY = 3;

  my $coefY    = ( $self->{H} - $self->{MARGIN} ) / $maxY; 

  
  my @col = ("blue", "red");
  my $inc = 1;
  my $i   = 0;
  
  foreach my $r (@{ $self->{PROFILES} }) {
    
    my $xp = 0;
    my $yp = 0;
    
    $self->{P}->setcolour($col[$i]);
    
    my @pol = ($self->{XBASE},$self->{YBASE});

    for (my $i=0; $i<$n; $i++) {
      
      my $x = $self->{XBASE} + ($r->[$i]->[0] - $minX)    * $coefX;
      my $y = $inc + $self->{YBASE} + Sets::max(0, $self->{ADD} + $r->[$i]->[1]) * $coefY ;
      #print "$x\t$y\n";
      #$self->{P}->line($xp,$yp,$x,$y) if ($i > 0);
      
      push @pol, ($x,$y);

      $xp = $x;
      $yp = $y;
      
      #print "$x\t$y\n";
    }
    
    push @pol, ($self->{XBASE}+$self->{W},$self->{YBASE});

    push @pol, $pol[0];
    push @pol, $pol[1];

    $self->{P}->polygon({ filled => 1}, @pol);


    #$inc ++;
    $i ++;
  }

  $self->{P}->setcolour("black");
  
  $self->{P}->line( $self->{XBASE}, $self->{YBASE},
		    $self->{XBASE} + $self->{W}, $self->{YBASE} );

  if (defined($self->{HLINE})) {

    $self->{P}->line( 
		     $self->{XBASE},
		     $self->{YBASE} + $self->{HLINE} * $coefY,
		     $self->{XBASE} + $self->{W}, 
		     $self->{YBASE} + $self->{HLINE} * $coefY
		    );
    
  }
  
  $self->{P}->line( $self->{XBASE}, $self->{YBASE},
		    $self->{XBASE}, $self->{YBASE} + $self->{H} );

  $self->{P}->setfont("Times", 8);
  
  my $k = 0;
  while ($k < $self->{H}/$coefY) {
    
    my $y = $k * $coefY;

    $self->{P}->line( $self->{XBASE},      $self->{YBASE} + $y,
		      $self->{XBASE} - 5 , $self->{YBASE} + $y);
    
    $self->{P}->text({ align => "right" }, $self->{XBASE} - 5, $self->{YBASE} + $y - 3, "$k ");
    $k ++;
  }

}

sub add {
  my ($self, $a) = @_;
  $self->{ADD} = $a;
}

sub output {

  my ($self, $f) = @_;

  if (defined($f)) {
    $self->{OUTFILE} = $f;
  }

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
