use lib "$ENV{FIREDIR}/SCRIPTS";
use lib "$ENV{FIREDIR}/SCRIPTS/PostScript-Simple-0.07/lib";

use PostScript::Simple;
use Table;
use Sets;
use Getopt::Long;
use Fasta;
use strict;

my $scriptdir        = "$ENV{FIREDIR}/SCRIPTS";

my $profiles         = undef;
my $summaryfile      = undef;
my $epsfile          = undef;
my $seqlen           = undef;
my $ps2pdf           = 0;
my $leftlabel        = 'leftlabel';
my $rightlabel       = 'rightlabel';
my $motifs_i         = undef;
my $motifs_k         = undef;
my $motifs_m         = undef;
my $clusters         = undef;
my $expfile          = undef;
my $outeps           = undef;
my $rna              = undef;
my $intersection     = 0;
my $overrep          = 0;
my $h                = 10;
my $fullmatrixfile   = undef;
my $forcea4          = 0;
my $fillup           = 0;
my $fastafile        = undef;

if (@ARGV == 0) {
  die "Usage: perl mi_draw_motif_map.pl --expfile=FILE --motifs=STR --clusters=STR --summaryfile=FILE --profiles=FILE --epsfile=FILE --seqlen=INT --ps2pdf=INT --outeps=FILE --rna=INT --motifs_k=STR --motifs_i=STR --motifs_m=STR --fastafile\n";
}

GetOptions ('summaryfile=s'   => \$summaryfile,
	    'fullmatrixfile=s'=> \$fullmatrixfile,
	    'profiles=s'      => \$profiles,
	    'epsfile=s'       => \$epsfile,
	    'expfile=s'       => \$expfile,
	    'rightlabel=s'    => \$rightlabel,
	    'leftlabel=s'     => \$leftlabel,
	    'seqlen=s'        => \$seqlen,
	    'overrep=s'       => \$overrep,
	    'forcea4=s'       => \$forcea4,
	    'motifs_i=s'      => \$motifs_i,
	    'motifs_m=s'      => \$motifs_m,
	    'motifs_k=s'      => \$motifs_k,
	    'h=s'             => \$h,
	    'clusters=s'      => \$clusters,
	    'fillup=s'        => \$fillup,
	    'outeps=s'        => \$outeps,
	    'intersection=s'  => \$intersection,
	    'rna=s'           => \$rna,
	    'fastafile=s'     => \$fastafile,
	    'ps2pdf=s'        => \$ps2pdf);


my $ta = Table->new;

#
# get the motifs, one way or another
#

my @a_motifs_i = undef;
if (defined($motifs_i)) {
  @a_motifs_i = split /\,/, $motifs_i;
}

my @a_motifs_k = undef;
if (defined($motifs_k)) {
  @a_motifs_k = split /\,/, $motifs_k;
}

my @a_motifs_m = undef;
if (defined($motifs_m)) {
  @a_motifs_m = split /\,/, $motifs_m;
}


my @a_clusters = split /\,/, $clusters;

#my @a_fillup = split /\,/, $fillup;


#
# load expfile
#
$ta->loadFile($expfile);
my $h_ref_exp = $ta->getIndexKV(0,1);

#
# diff lengths
#

my $h_ref_len  = undef;
my $max_seqlen = -1;
if (defined($fastafile) && !defined($seqlen)) {

  $h_ref_len = {};
  my $fa = Fasta->new;
  $fa->setFile($fastafile);
  while (my $a_ref = $fa->nextSeq()) {
    my ($n, $s) = @$a_ref;
    my $l = length($s);
    $h_ref_len->{$n} = $l;
    if ($l > $max_seqlen) {
      $max_seqlen = $l;
    }
  }

  if ($rna == 1) {
    $leftlabel  = 0;
    my $tmp = $max_seqlen; $tmp --;
    $rightlabel = "+$tmp";
  }

} elsif (defined($seqlen)) {
  $max_seqlen = $seqlen;

} else {
  die "Exiting .. please define seqlen.\n";

}





#
# load summary file
#

my %h_sum     = ();
my %h_col     = ();
my @colors = ('red', 'blue', 'green');
my $cntm      = 0;

my %overclu    = ();
my $a_ref_sum  = undef;
my $h_ref_sum8 = undef;
my $h_ref_sum0 = undef;
if (defined($summaryfile)) {
  $ta->loadFile($summaryfile);
  $a_ref_sum  = $ta->getArray();
  $h_ref_sum8 = $ta->getIndex(8);
  $h_ref_sum0 = $ta->getIndex(0);
  
}


my %pvalues = ();
if (defined($fullmatrixfile)) {
  #
  # load pvalue matrix file
  #
  $ta->loadFile($fullmatrixfile);
  my $a_ref = $ta->getArray();
  shift @$a_ref;

  foreach my $r (@$a_ref) {
    my $m = shift @$r;
    $pvalues{ $m } = $r;
  }

}



if (!defined($motifs_m) && (defined($summaryfile))) {

  if (defined($motifs_i)) {
    
    # i
    foreach my $i (@a_motifs_i) {
      my $r = $a_ref_sum->[$i];
      $h_sum{ $r->[0] } = $r;
      $h_col{ $r->[0] } = $colors[$cntm++];
    }

  } elsif (defined($motifs_k)) {
    
    # kmer
    foreach my $k (@a_motifs_k) {
      my $r = $h_ref_sum8->{$k};
      $h_sum{ $r->[0] } = $r;
      $h_col{ $r->[0] } = $colors[$cntm++];
    }
  }

}

my $prefix = int(rand(10000));

if (defined($motifs_m)) {

  # motif
  die "Please define -rna INT\n" if (!defined($rna));

  my $cnt = 0;
  foreach my $m (@a_motifs_m) {
    my $r = [$m, $rna];
    $h_sum{ $r->[0] } = $r;
    $h_col{ $r->[0] } = $colors[$cntm++];

    if ($overrep == 1) {
      my $nn = @{ $h_ref_sum0->{$m} };
      for (my $i=12; $i<$nn; $i++) {
	$overclu{ $h_ref_sum0->{$m}->[$i] } = 1;
      }
    }
    my $mo = Sets::myre2wm($m);
    if ($rna == 1) {
      $mo =~ s/T/U/g;
    }

    my $tmptxt = "$prefix$cnt.txt";
    my $tmpeps = "$prefix$cnt.eps";
    open OUT, ">$tmptxt" or die "cannot open $tmptxt\n";
    print OUT $mo;
    close OUT;
    print "Outputing $m\n";
    system("$scriptdir/weblogo/seqlogo -f $tmptxt -F EPS  -a -c -M -n -Y -w 5 -h 3 > $tmpeps");
    $cnt ++;
  }
}

#
# load clusters
#
my %h_clu = ();
foreach my $r (@a_clusters) {
  $h_clu{ $r } = $r;
}


#
# load profiles, record gene => 
#
my %CLU_GENE_MOT = ();
my %H            = ();

$ta->loadFile($profiles);
my $a_ref = $ta->getArray();
foreach my $r (@$a_ref) {
  
  #print join("\t", @$r) . "\n";
  
  # cluster
  my $c = $h_clu{ $h_ref_exp->{ $r->[1] } };
  #print "Cluster $c ($h_ref_exp->{$r->[1]})  ($r->[1])\n";

  # gene must be in data
  next if (!defined($h_ref_exp->{ $r->[1] }));

  # must be a good cluster
  next unless (($clusters eq "-1") || defined( $c ));

  # must be a good motif
  next unless (defined($h_sum{$r->[0]}));

  my @a = ($r->[2], $r->[3]);
  push@{ $CLU_GENE_MOT{$h_ref_exp->{ $r->[1] }}{$r->[1]}{$r->[0]} }, \@a;

  $H{ $r->[1] } = 1;


}

#
# how many genes do we have
#
my $n = 0;

my @a_sorted_clusters = ();
if (defined($fullmatrixfile)) {

  my $nc = scalar( @{$pvalues{$a_motifs_m[0]}} );
  
  my @scores = ();
  
  # go thru all clusters
  for (my $i=0; $i<$nc; $i++) {
    my $minp = 100000;
    my $minm = undef;

    # go thru all motifs
    foreach my $m (@a_motifs_m) {
      #print "   $pvalues{$m}->[$i]\n";
      if ($pvalues{$m}->[$i] < $minp) {
	$minp = $pvalues{$m}->[$i];
	$minm = $m;
      }
    }
    $scores[$i] = $minp;

    print "$i\t$minp\t$minm\n";

  }


  #foreach my $m (@a_motifs_m) {
  # my $a_ref_p = $pvalues{$m};
  #  for (my $i=0; $i<@$a_ref_p; $i++) {
  #    $scores[$i] += $a_ref_p->[$i];
  #  }
  #}

  my $a_ref = Sets::order(\@scores, 0);
  @a_sorted_clusters = @$a_ref;

  my $nn = @a_sorted_clusters;

  if ($fillup>0) {
    
    my @newa = ();
   
    for (my $i=0; $i<$fillup; $i++) {
      push @newa, shift(@a_sorted_clusters);
    }	
    
    push @newa, reverse(@a_sorted_clusters);
    
    @a_sorted_clusters = @newa;

    #<STDIN>;
  }
    
  


} else {

  @a_sorted_clusters = sort { $a <=> $b } (keys(%CLU_GENE_MOT));

}

foreach my $c (@a_sorted_clusters) {
    
  if (defined($overrep) && ($overrep == 1) && (!defined($overclu{$c}))) {
    next;
  }
  
  foreach my $g (sort(keys(%{$CLU_GENE_MOT{$c}}))) {
    if ($intersection == 1) {
      my @tmp = keys(%{$CLU_GENE_MOT{$c}{$g}});
      my $nm = @tmp;
      if ($nm != $cntm) {
	next;
      }
    }
    $n ++;
  }
}


#die "Not enough genes.\n" if ($n < 50);

my $ybase  = 70;
my $xright = 50;
my $xleft  = 30;
my $xsize  = 595;
my $xscale = ($xsize - $xright - $xleft) / $max_seqlen;


my $ysize = $h * $n + $ybase + 100;
if ($forcea4 == 1) {
  $ysize = Sets::min(842, $ysize);
}

my $p = new PostScript::Simple(xsize     => $xsize,
			       ysize     => $ysize,
			       
			       colour    => 1,
			       eps       => 1,
			       units     => "pt");

$p->setlinewidth(0.5);



if (defined($epsfile)) {
  my $e  = new PostScript::Simple::EPS(file => $epsfile);
  $p->_add_eps($e, $xsize/2 - 80,  $ysize - $ybase - 20); # add eps object to postscript object again
}

my $pos = $xsize / 2 - $cntm * 80;
for (my $i=0; $i<$cntm; $i++) {
  my $tmpeps = "$prefix$i.eps";
  my $e  = new PostScript::Simple::EPS(file => $tmpeps);
  die "Cannot find $tmpeps\n" if (! -e $tmpeps);
  $p->_add_eps($e, $pos + $i * 120,  $ysize - $ybase - 20); # add eps object to postscript object again
  
  draw_symbol($pos + $i * 120 + 20,  $ysize - 10, 'otriangle', 1, $colors[$i], Sets::max(2,int($h/2)));

}

$p->setcolour("black");
$p->setfont("Courrier", 8);

$p->text({ align => 'left' }, $xsize - ($xright - 1), $ysize - ($ybase  - 8), $rightlabel);
$p->text({ align => 'center' }, $xleft, $ysize - ($ybase  - 8), $leftlabel );

my $i = 0;

#
# go through all clusters
#

my $cnt_c = 0;
foreach my $c (@a_sorted_clusters) {

  #
  # go through all genes
  #

  if (defined($overrep) && ($overrep == 1) && (!defined($overclu{$c}))) {
    next;
  }
  
  print "Cluster $c.\n";


  my $nbgenes = 0;
  foreach my $g (sort(keys(%{$CLU_GENE_MOT{$c}}))) {
    if ($intersection == 1) {
      my @tmp = keys(%{$CLU_GENE_MOT{$c}{$g}});
      my $nm = @tmp;
      if ($nm != $cntm) {
	next;
      }
    }
    $nbgenes ++;
  }

  if ($nbgenes == 0) {
    next;
  }
  
  if ( ($ybase  + ($i+$nbgenes) * $h) > $ysize) {
    last;
  }

  $p->setcolour("black");
  $p->setfont("Courrier", 8);  
  $p->text({ align => 'left' }, 2, $ysize - ($ybase  + ( $i + ($nbgenes / 2) ) * $h), "C$c");

  foreach my $g (sort(keys(%{$CLU_GENE_MOT{$c}}))) {

    # may draw only if the total number of motif is right
    if ($intersection == 1) {

      my @tmp = keys(%{$CLU_GENE_MOT{$c}{$g}});
      my $nm = @tmp;

      if ($nm != $cntm) {
	next;
      }

	
    }
    

    print "  Gene $g.\n";

    $p->setcolour("black"); $p->setfont("Courrier", 6);
    $p->text({ align => 'left' }, $xsize - ($xright - 5), $ysize - ($ybase  + $i * $h + $h/2), $g);
    
    my $x1 = undef;
    my $x2 = undef;
    if ($rna == 0) {

      $x1 = $xsize -   $xright;
      my $this_seqlen = undef;
      if (!defined($h_ref_len)) {
	$this_seqlen = $seqlen;
      } else {
	$this_seqlen = $h_ref_len->{$g};
      }

      $x2 = $xsize - ( $xright + $xscale * $this_seqlen );

    } else {
      
      $x1 = $xleft; 
      my $this_seqlen = undef;
      if (!defined($h_ref_len)) {
	$this_seqlen = $seqlen;
      } else {
	$this_seqlen = $h_ref_len->{$g};
      }

      $x2 = $xleft + $this_seqlen * $xscale; 

      
    }
    
    $p->line($x1, $ysize - ($ybase + $i * $h),
	     $x2, $ysize - ($ybase + $i * $h)  );

    #
    # go through all motifs
    #

    foreach my $m (keys(%{$CLU_GENE_MOT{$c}{$g}})) {

      print "    Motif $m.\n";

      if (!defined($rna)) {
	print "Getting RNA/DNA info from summaryfile.\n";
	$rna = $h_sum{$m}->[1];
      }
	
      my $col = $h_col{$m};

      

      foreach my $r (@{$CLU_GENE_MOT{$c}{$g}{$m}}) {
	
	my ($p, $s) = @$r;
	my $y = $ysize - ($ybase + $i * $h);
	my $x = undef;
	
	#print "s=$s\n";
	
	if ($rna == 0) {
	  $x = $xsize - ($xright + $p * $xscale);
	} else {
	  $x = $xleft + $p * $xscale;
	}

	if ($s == -1) {
	  #$col = "black";
	}
	draw_symbol($x, $y, 'otriangle', $s, $col, Sets::max(2,int($h/2)));

      }

    }

    $i++;
  }
  
  $cnt_c ++;
  
  if (($fillup > 0) && ($fillup == $cnt_c)) {
    
    $p->setlinewidth(2);
    $p->setcolour("black");
    $p->line(0+20, $ysize - ($ybase  + $i * $h + 2), $xsize-15 , $ysize - ($ybase + $i * $h + 2));
    
    $p->setlinewidth(0.5);

    
  }

  
  
  $i+=2;

 
}


sub draw_symbol {
  my ($x,$y,$t,$s,$c,$w) = @_;

  print "Draw $c at $x,$y,$s,$w\n";

  $p->setcolour($c);

  if (!defined($w)) {
    $w = 3;
  }

  #if ($s == 1) {
  #  $p->setcolour("blue");
  #} else {

  #}

  my @pol = undef;
  
  
  

  if ($t eq 'triangle') {
    @pol = ($x, $y-$w,
	    $x-$w, $y+$w,
	    $x+$w, $y+$w);

  } elsif ($t eq 'otriangle') {
    
    if ($s == 1) {
      @pol = ( $x-$w, $y + $w,
	       $x-$w, $y - $w,
	       $x+$w, $y, 
	       $x-$w, $y + $w
	     );
    } else {
      @pol = ( $x+$w, $y + $w,
	       $x+$w, $y - $w,
	       $x-$w, $y,
	       $x+$w, $y + $w);
    }
  } elsif ($t eq 'square') {

    @pol = ( $x - 2, $y + $w,
	     $x + 2, $y + $w,
	     $x + 2, $y - $w,
	     $x - 2, $y - $w);
  }
    
  
  
  $p->polygon({ filled => 1}, @pol);
  
}

print "Producing $outeps\n";
$p->output("$outeps");


if ($ps2pdf == 1) {
  my $outpdf = $outeps; $outpdf =~ s/\.eps/\.pdf/;
  system("ps2pdf -dEPSCrop -dAutoRotatePages=/None $outeps $outpdf");
}

for (my $i=0; $i<$cntm; $i++) {
  my $tmpeps = "$prefix$i.eps";
  my $tmptxt = "$prefix$i.txt";
  unlink $tmpeps; unlink $tmptxt;
}
