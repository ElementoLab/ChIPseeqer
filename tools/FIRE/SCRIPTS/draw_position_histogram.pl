#
# drawing motif position histograms
#
# Most of this script was written by Danny Lieber !
#


use lib "$ENV{FIREDIR}/SCRIPTS";
use lib "$ENV{FIREDIR}/SCRIPTS/PostScript-Simple-0.07/lib";

use Getopt::Long;
use Table;
use Sets;
use Data::Dumper ;
use PostScript::Simple;

use strict;

my $nbclusters = 0 ;
my $expfile     = undef;
my $data        = undef;
my $summaryfile = undef;
my $pvaluematrixfile = undef ;
my $ps2pdf      = 1;

my $bins        = 10;


my $targetfile = undef;
my $profiles = undef;
my $expfile = undef;

my $outeps = "clusters.eps";
my $outpdf = "clusters.pdf";
my $outpng = "clusters.png";
my $outdir = "./";

my $seqlen = undef;
my $rna    = undef;
my $onlybins = undef;

if (@ARGV == 0) {
    die "Usage: perl draw_position_histogram.pl --profiles=FILE --targetfile=FILE --outdir=DIR\n";
}
GetOptions ( 'expfile=s'            => \$expfile,
'profiles=s'           => \$profiles,
'targetfile=s'         => \$targetfile,
'outdir=s'             => \$outdir,
'rna=s'                => \$rna,
'seqlen=s'             => \$seqlen,
'onlybins=s'           => \$onlybins,
'bins=s'               => \$bins
);

if (! -e $outdir) { 
	print "Outdir $outdir does not exist. Creating...";
	system("mkdir -p $outdir\n");
	print "Done\n";
}

#
# 
#
my @bins_set = ();
if (defined($onlybins)) {
	@bins_set = split /\,/, $onlybins;
}

my @MOTIFS       = ();
my %TARGETS      = ();
my %POS_ENRICHED = ();
my %POS_OTHER    = ();

my $ta = Table->new;
$ta->loadFile($targetfile);
my $a_ref_sum = $ta->getArray();
my $nbcols = $ta->getNbColumns;
my $nbrows = $ta->getNbRows;

foreach my $r (@$a_ref_sum) {
	push @MOTIFS, $r->[0];
	for (my $i=1; $i<@$r; $i++) {
		$TARGETS{$r->[0]}{$r->[$i]} = 1; # hash motif-gene
	}
}


#
#  read in expfile
#
$ta->loadFile($expfile);
my $h_ref_exp = $ta->getIndexKV(0,1);
my $a_ref_exp = $ta->getColumn(0);



### load profiles
my $ta = Table->new;       
$ta->loadFile($profiles);
my $a_ref_pro = $ta->getArray();

foreach my $r (@$a_ref_pro) {
	
	# get cluster
	my $c = $h_ref_exp->{ $r->[1] };
	next if (!defined($c));
	
	# get relative position
	my $percentage = undef;
	if ($rna == 0) {
		$percentage = sprintf("%2.2f", ($seqlen - $r->[2]) / $seqlen);
	} else {
		$percentage = sprintf("%2.2f", $r->[2] / $seqlen);
	}
	
	# add it to the right set depending on whether it is a target or not
	if(exists $TARGETS{$r->[0]}{$r->[1]}){
		if (!defined($onlybins) || (defined($onlybins) && Sets::in_array($c, @bins_set))) { 
			push @{$POS_ENRICHED{$r->[0]}}, $percentage;  
		}
		
	} else {
		push @{$POS_OTHER{$r->[0]}}, $percentage;
	}
	
}

#
# now iterates thru the motifs, and draw an histogram for each of them
#
foreach my $motif (@MOTIFS){
	
	next if (!defined @{$POS_ENRICHED{$motif}});
	
	print $motif, "\n";
	
	my $enriched_hist_ref = _histogram($POS_ENRICHED{$motif}, $bins);
	my $other_hist_ref    = _histogram($POS_OTHER{$motif},    $bins);
	
	next if (!defined($other_hist_ref));
	
	my @enriched_hist     = @{$enriched_hist_ref};
	my @other_hist        = @{$other_hist_ref};
    
	my $max_enriched      = 0;
	my $sum_enriched      = 0;
	
	my $max_other         = 0;
	my $sum_other         = 0;
	
	my @enriched_dist     = normalize($enriched_hist_ref,\$sum_enriched, \$max_enriched);
	my @other_dist        = normalize($other_hist_ref,\$sum_other, \$max_other);
	
	print "Targets\tSum: $sum_enriched\tMax: $max_enriched\n";
	print "Non-targets\tSum: $sum_other   \tMax: $max_other\n";
	
	foreach my $p (@enriched_dist){
		print sprintf("%2.2f",$p);
		print "\t";
	}
	print "\n";
	foreach my $p (@other_dist){
		print sprintf("%2.2f",$p);
		print "\t"; 
	}
	print "\n";
	
	#
	#  START DRAWING
	#
	#
	
	my $xsize = 600 ;
	my $ysize = 600;
	
	my $xbase = 100;
	my $ybase = 50;
	
	my $barwidth = ($xsize-2*$xbase)/$bins;
	$barwidth = min2( 50, $barwidth);
	my $w  = $barwidth+2;
	
	my $height = 200;
	my $scale = max2 ($max_enriched/$sum_enriched, $max_other/$sum_other);
	
	#my $maxheight = $height * $scale ;
	
	
	
	my $xoffset = 50;
	my $yoffset = $height * $max_other/ $sum_other / $scale + 30;
	
	
	my $p          = new PostScript::Simple(xsize     => $xsize,
	ysize     => $ysize,
	colour    => 1,
	eps       => 1,
	units     => "pt");
	
	$p->setlinewidth(1.0);
	
	# Write Titles
	$p->setcolour("red");
	$p->setfont("Garamond", 30);
	$p->text({ align => "center" }, $xbase+$xoffset+$bins*$w/2, $ybase+$yoffset+$height*$max_enriched/ $sum_enriched / $scale+22,  $motif);
	
	$p->setcolour("blue");
	$p->setfont("Garamond", 12);
	$p->text({ align => "center" }, $xbase+$xoffset+$bins*$w/2+20, $ybase+$yoffset+$height*$max_enriched/ $sum_enriched / $scale+5,"Targets: $sum_enriched\nNon-targets: $sum_other");
	
	$p->setcolour("black");
	$p->setfont("Garamond", 12);
	$p->text({ align => "center" }, $xbase+$xoffset+$bins*$w/2, $ybase-30,"Position (\% of sequence length)");
	$p->text({rotate => 90}, $xbase+10, $ybase+.75*$height, "Frequency");
	
	my $border = 3;
	
	# draw rest of histogram
	
	# enriched
	for (my $i=0; $i<$bins; $i++){
		$p->setcolour("blue");
		#  print "$enriched_dist[$i]\t";
		$p->box({filled => 1},$xbase+$xoffset+$i*$w, $ybase+$yoffset, $xbase+$xoffset+$i*$w+$barwidth,$ybase+$yoffset+$height*$enriched_dist[$i]/$scale);
		# top left shadow
		$p->setcolour(152,192,255);  #light blue
		$p->polygon({filled => 1},
		$xbase+$xoffset+$i*$w, $ybase+$yoffset,
		$xbase+$xoffset+$i*$w+$border,$ybase+$yoffset+$border,
		$xbase+$xoffset+$i*$w+$border,$ybase+$yoffset+$height*$enriched_dist[$i]/$scale-$border,
		$xbase+$xoffset+$i*$w+$barwidth-$border,$ybase+$yoffset+$height*$enriched_dist[$i]/$scale-$border,
		$xbase+$xoffset+$i*$w+$barwidth,$ybase+$yoffset+$height*$enriched_dist[$i]/$scale,
		$xbase+$xoffset+$i*$w,$ybase+$yoffset+$height*$enriched_dist[$i]/$scale,
		);
		
		#bottom right shadow
		$p->setcolour("darkblue");  #dark blue
		$p->polygon({filled => 1},
		$xbase+$xoffset+$i*$w, $ybase+$yoffset,
		$xbase+$xoffset+$i*$w+$border,$ybase+$yoffset+$border,
		$xbase+$xoffset+$i*$w+$barwidth-$border,$ybase+$yoffset+$border,
		$xbase+$xoffset+$i*$w+$barwidth-$border,$ybase+$yoffset+$height*$enriched_dist[$i]/$scale-$border,
		$xbase+$xoffset+$i*$w+$barwidth,$ybase+$yoffset+$height*$enriched_dist[$i]/$scale,
		$xbase+$xoffset+$i*$w+$barwidth,$ybase+$yoffset,
		);
		
		$p->setcolour("blue");
		$p->setfont("Garamond", 12);
		if ($i == 0){
			$p->text({ align => "center" }, $xbase+$xoffset,  $ybase+$yoffset-10, 0);
		}
		$p->text({ align => "center" }, $xbase+$i*$w+$xoffset+$barwidth,  $ybase+$yoffset-10, sprintf("%0.2f",($i+1)/$bins)+0);
		$p->setcolour("red");  $p->setfont("Arial", 18);
		$p->text({ align => "right" }, $xbase+$xoffset-30,  $ybase+$yoffset+$height/2.2, "Targets")
    }
	
	# other
	for (my $i=0; $i<$bins; $i++){
		
		
		$p->setcolour("blue");
		#  print "$other_dist[$i]\t";
		$p->box({filled => 1},$xbase+$xoffset+$i*$w, $ybase, $xbase+$xoffset+$i*$w+$barwidth,$ybase+$height*$other_dist[$i]/$scale);
		
		# top left shadow
		$p->setcolour(152,192,255);  #light blue
		$p->polygon({filled => 1},
		$xbase+$xoffset+$i*$w, $ybase,
		$xbase+$xoffset+$i*$w+$border,$ybase+$border,
		$xbase+$xoffset+$i*$w+$border,$ybase+$height*$other_dist[$i]/$scale-$border,
		$xbase+$xoffset+$i*$w+$barwidth-$border,$ybase+$height*$other_dist[$i]/$scale-$border,
		$xbase+$xoffset+$i*$w+$barwidth,$ybase+$height*$other_dist[$i]/$scale,
		$xbase+$xoffset+$i*$w,$ybase+$height*$other_dist[$i]/$scale,
		);
		
		#bottom right shadow
		$p->setcolour("darkblue");  #dark blue
		$p->polygon({filled => 1},
		$xbase+$xoffset+$i*$w, $ybase,
		$xbase+$xoffset+$i*$w+$border,$ybase+$border,
		$xbase+$xoffset+$i*$w+$barwidth-$border,$ybase+$border,
		$xbase+$xoffset+$i*$w+$barwidth-$border,$ybase+$height*$other_dist[$i]/$scale-$border,
		$xbase+$xoffset+$i*$w+$barwidth,$ybase+$height*$other_dist[$i]/$scale,
		$xbase+$xoffset+$i*$w+$barwidth,$ybase,
		);
		
		$p->setfont("Garamond", 12);
		if ($i == 0){
			$p->text({ align => "center" }, $xbase+$xoffset,  $ybase-10, 0);
		}
		$p->text({ align => "center" }, $xbase+$i*$w+$xoffset+$barwidth,  $ybase-10, sprintf("%0.2f",($i+1)/$bins)+0);
		$p->setcolour("red");  $p->setfont("Arial", 18);
		$p->text({ align => "right" }, $xbase+$xoffset-30,  $ybase+$height/2.2, "Non-targets")
	}
	
	## draw tick marks
	$p->setcolour("black");
	$p->setlinewidth(1);
	$p->setfont("Garamond", 10);
	
	$p->line($xbase+$xoffset, $ybase+$yoffset, $xbase + $xoffset, $ybase+$yoffset+$height * $max_enriched / $sum_enriched  / $scale);
	$p->line($xbase+$xoffset, $ybase, $xbase + $xoffset, $ybase+ $height * $max_other / $sum_other  / $scale);
	
	
	my $i = 0;
	while ($i <= $height * $max_enriched / $sum_enriched  / $scale){
		#  print $i, "\n";
		$p->line($xbase+$xoffset-5, $ybase+$yoffset+$i+.2, $xbase + $xoffset+ 1, $ybase+$yoffset+$i + .2);
		$p->text({ align => "right" }, $xbase+$xoffset-6,  $ybase+$yoffset+$i, sprintf("%0.2f", $i*$scale/$height)+0);
		$i += .05 * $height/$scale;
		#   print "Next i: $i\n";
	}
	my $i = 0;
	while ($i <= $height * $max_other / $sum_other  / $scale){
        #    print $i, "\n";
		$p->line($xbase+$xoffset-5, $ybase+$i+.2, $xbase + $xoffset+ 1, $ybase+$i + .2);
		$p->text({ align => "right" }, $xbase+$xoffset-6,  $ybase+$i, sprintf("%0.2f", $i*$scale/$height)+0);
		$i += .05 * $height/$scale;
		#    print "Next i: $i\n";
	}
	
	$outeps = "$outdir/".$motif."_histogram.eps";
	$outpdf = $outeps;
	$outpng  = $outeps;
	$outpdf =~ s/eps/pdf/;  
    $outpng =~ s/eps/png/;  
	
	$p->output($outeps);
	
	if ($ps2pdf == 1) {
		print "Creating PDF $outpdf  ... ";
		system("ps2pdf -dEPSCrop -dAutoRotatePages=/None $outeps $outpdf");
		#system("pstoimg -crop a -antialias $outeps ");
		print "Done.\n";
	}   
	
	
}

sub normalize {
    my($ref, $sum, $max) = @_;
    my @data = @$ref;
    my @norm = ();
    $$sum = 0;
    $$max = 0;
    foreach (@data){
        $$sum += $_; 
        $$max = ( $_ > $$max )? $_ : $$max;
    }
    #  print "Sum $$sum \n";
    foreach (@data){
        push @norm, $_ / $$sum; 
    }
    
    return @norm;
}

sub max2 ($$) { $_[$_[0] < $_[1]] }
sub min2 ($$) { $_[$_[0] > $_[1]] }


#------------------------------------------------------------#
######################################
#    Courtesy of Snehanshu Shah      #
#        GD::Graph::histogram        #
######################################

# returns reference to histogram array
sub _histogram { 
	my ($data, $nbins ) = @_;
	
	my $cp = _bins( $data, $nbins);
	my $binArrRef = _histogram_frequency( $data, $cp );
	
	
	return $binArrRef;
}

###################################################################
# _bins - calculates the bins usings Scott's algorithm
#
# Arguements:
#
# $data  (Vector)  - Data values
#
# $nbins (Integer) - Number of bins to create. If $nbins is undef
#                    the number of bins is calculated using Scott's
#                    algorithm
#
###################################################################
sub _bins {
	my ( $data, $nbins ) = @_;
	
	if( !defined $data ) { return; }
	
	my $calcBins = ( defined $nbins )? 0 : 1;
	my $cnt = 0;
	my $mean= 0;
	my $max = my $min = $data->[0];
	foreach (@$data) {
		$mean += $_;
		$min = ( $_ < $min )? $_ : $min;
		$max = ( $_ > $max )? $_ : $max;
		$cnt++;
	}
	$mean /= $cnt if( $cnt > 1 );
	
	my $sumsq = 0;
	$nbins = 1 if( $calcBins );
	my $s = 0;
	if( $cnt > 1 ) {
		foreach (@$data) {
			$sumsq += ( $_ - $mean )**2;
		}
		$s = sqrt( $sumsq / ($cnt - 1));
		$nbins = 3.49 * $s / $cnt**0.33 if( $s > 0 && $calcBins );
	}
	
	my $binwidth = ( $max - $min ) / $nbins;
	
	my $lower = $min;
	my $upper = $lower;
	
	my $bins;
	my @cutPoints;
	my $cntr = 0;
	while ( $upper <= $max && $cntr < $nbins) {
		$upper = $lower + $binwidth;
		push( @cutPoints, [$lower, $upper] );
		$lower = $upper;
		$cntr++;
	}
	
	return \@cutPoints;
}

###################################################################
# _histogram_frequency - bins the data
#
#     Lower Boundry <= data value < Upper Boundry
#
# Arguements:
#
# $data  (Vector)  - Data values
#
# $nbins (Integer) - Vector containing the cutpoints to bin the data
#
###################################################################
sub _histogram_frequency {
	my ( $data, $cutPoints ) = @_;
	
	if( !defined $data || !defined $cutPoints ) { return; }
    
	my @freqs;
	foreach (@$cutPoints) {
		push( @freqs, 0 );
	}
	
	foreach (@$data) 
    {
		for( my $i = 0; $i < scalar( @$cutPoints ); $i++ ) 
		{
			if( ($_ >= $cutPoints->[$i]->[0] && $_ < $cutPoints->[$i]->[1])
				||
				($i == (scalar (@$cutPoints) - 1) && $_ >= $cutPoints->[$i]->[1]) ) 
			{	
				
				$freqs[$i]++;
			}
		}
    }
	return \@freqs;
}