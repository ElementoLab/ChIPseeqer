package Subs ;

use lib "$ENV{FIREDIR}/SCRIPTS";
use lib "$ENV{FIREDIR}/SCRIPTS/PostScript-Simple-0.07/lib";

use Algorithm::Cluster ;
use PostScript::Simple ;
use Table ;
use Sets ;
use Data::Dumper ;
use strict ;

sub varnormTable
{
    my $ref_M =  shift @_ ;
    my $a_ta_M = copyTable($ref_M) ;
    my $a_ta_H = shift(@$a_ta_M) ;
    my $dummy = shift(@$a_ta_H) ;
    my $a_ta_R ;
    foreach my $r (@$a_ta_M)
    {
        push(@$a_ta_R, $r->[0]) ;
        shift @$r ;
    }
    foreach my $r (@$a_ta_M)
    {
	my $mean = Sets::average($r) ;
	my $std = Sets::stddev($r) ;
	for (my $i=0 ; $i<@$a_ta_H ; $i++)
	{
	    $r->[$i] = ($r->[$i]-$mean)/$std ;
	}
    }
    for (my $i=0 ; $i<@$a_ta_R ; $i++)
    {
        unshift(@{$a_ta_M->[$i]}, $a_ta_R->[$i]) ;
    }
    unshift (@$a_ta_H, $dummy) ;
    unshift (@$a_ta_M, $a_ta_H) ;
    return $a_ta_M ;
}

sub copyTable
{
    my $a_ta_M = shift @_ ;
    my $new_M ;
    for (my $i=0 ; $i<@$a_ta_M ; $i++)
    {
	my $r = $a_ta_M->[$i] ;
	for (my $j=0 ; $j<@$r ; $j++)
	{
	    $new_M->[$i]->[$j] = $r->[$j] ;
	}
    }
    return $new_M ;
}

sub transposeTable
{
    my $ref_M =  shift @_ ;
    my $a_ta_M = copyTable($ref_M) ;
    my $a_ta_H = shift(@$a_ta_M) ;
    my $dummy = shift(@$a_ta_H) ;
    my $a_ta_R ;
    foreach my $r (@$a_ta_M)
    {
        push(@$a_ta_R, $r->[0]) ;
        shift @$r ;
    }

    my $new_M ;
    for (my $i=0 ; $i<@$a_ta_H ; $i++)
    {
	for (my $j=0 ; $j<@$a_ta_R ; $j++)
	{
	    $new_M->[$i]->[$j] = $a_ta_M->[$j]->[$i] ;
	}
    }
    for (my $i=0 ; $i<@$new_M ; $i++)
    {
	my $r = $new_M->[$i] ;
	unshift(@$r, $a_ta_H->[$i]) ;
    }
    unshift(@$a_ta_R, $dummy) ;
    unshift(@$new_M, $a_ta_R) ;
    
    return $new_M ;
}

sub printHeatmapBars
{
    my ($a_ta_M, $fn, $vline, $hline, $vcol, $min, $mid, $max, $h, $w, $xplus, $yplus, $header_motif, $row_motif, $s_y, $s_h, $s_w, $sep, $res, $scalefont, $uppertext, $lowertext,$lcol, $hcol, $mcol) = @_ ;

    my $a_ta_H = shift(@$a_ta_M) ;
    my $dummy = shift(@$a_ta_H) ;

    my $xbase        = 100 ;
    my $ybase        = 150 ;
    my $ysize        = $ybase + $h * scalar(@$a_ta_M) + $yplus ;
    my $xsize        = $xbase + $w * scalar(@$a_ta_H) + $xplus ;

    my $p = new PostScript::Simple(xsize     => $xsize,
                                   ysize     => $ysize,
                                   colour    => 1,
                                   eps       => 1,
                                   units     => "pt") ;
    $p->setlinewidth(0.5) ;

    if ($header_motif eq "true" or $row_motif eq "true")
    {
        system("mkdir ./Temp") if !(-d "./Temp") ;
    }

    for (my $j=0; $j<@$a_ta_H; $j++)
    {
        if ($header_motif eq "true")
        {
            my $motif = $a_ta_H->[$j] ;
	    $motif =~ /(\S+)\/(\S+)/ ;
	    $motif = $1 if defined $1 ;
	    my $name = $2 ;
	    $motif =~ s/\/$// ;
	    print $motif, "\t", $name, "\n" ;
            my $mo = Sets::myre2wm($motif);
            open OUT, "> ./Temp/$motif.txt" or die "cannot open\n";
            print OUT $mo;
            close OUT;
            system("../FIRE/SCRIPTS/weblogo/seqlogo -f ./Temp/$motif.txt -F EPS  -a -c -M -n -Y -w 5 -h 3 > ./Temp/$motif.eps");

            my $e  = new PostScript::Simple::EPS(file => "./Temp/$motif.eps");
            my $eh = $e->height;
            my $ew = $e->width;

            # height must be $h, so scale down to $h = k * LO
            $e->scale($w / $eh);
            $e->rotate(90);
            my $ew_new = int(0.5 + $ew * $h / $eh);
            $p->_add_eps($e, $xbase + ($j+1.1) * $w,  $ysize - ($ybase-3));
	    if (defined $name or $name ne "")
	    {
		$p->setcolour("black");
		$p->setfont("TimesBold", 10);
		$p->text( { rotate => 90 , align => "left"}, $xbase + $j * $w + 3*$w/4, $ysize - ($ybase-5-$w/$eh*$ew), $name);
	    }
	    else
	    {
		$p->setcolour("black");
                $p->setfont("TimesBold", 10);
                $p->text( { rotate => 90 , align => "left"}, $xbase + $j * $w + 3*$w/4, $ysize - ($ybase-5-$w/$eh*$ew), "-");
	    }
        }
        else
        {
            $p->setcolour("black");
            $p->setfont("Courrier", 6);
            $p->text( { rotate => 90 }, $xbase + $j * $w + 3*$w/4, $ysize - ($ybase-3), $a_ta_H->[$j]);
        }
    }

    my @col = () ;
    for (my $i=0; $i<@$a_ta_M; $i++)
    {
        my $r  = $a_ta_M->[$i] ;
        my $name = shift @$r;
        print $name,"\n" ;
        for (my $j=0; $j<@$r; $j++)
        {
            my $v = $r->[$j] ;

            #defining the color
            if (defined $mcol)
            {
                if ($v>$mid)
                {
                    @col = Sets::interp_general( $v, $mcol, $hcol, $mid, $max);
                }
                else
                {
                    @col = Sets::interp_general( $v, $lcol, $mcol, $min, $mid);
                }
            }
            else
            {
                @col = Sets::interp_general( $v, $lcol, $hcol, $min, $max);
            }

            $p->setcolour(@col);
            $p->box({filled => 1},
                    $xbase + $j * $w,      $ysize - ($ybase + $i*$h) ,
                    $xbase + $j * $w + $w, $ysize - ($ybase + ($i*$h+$h)));
	    
	    $p->setcolour(@$vcol);    
	    $p->line($xbase + $j * $w, $ysize - $ybase, $xbase + $j * $w, $ysize - ($ybase + ($i+1) * $h)) if $vline eq "true";
            $p->line($xbase , $ysize - ($ybase + $i*$h), $xbase + ($j+1) * $w, $ysize - ($ybase+$i*$h)) if $hline eq "true";
        }

        if ($row_motif eq "true")
        {
            my $motif = $name ;
            my $mo = Sets::myre2wm($motif);
            open OUT, "> ./Temp/$motif.txt" or die "cannot open\n";
            print OUT $mo;
            close OUT;
            system("../FIRE/SCRIPTS/weblogo/seqlogo -f ./Temp/$motif.txt -F EPS  -a -c -M -n -Y -w 5 -h 3 > ./Temp/$motif.eps");

            my $e  = new PostScript::Simple::EPS(file => "./Temp/$motif.eps");
            my $eh = $e->height;
            my $ew = $e->width;

            # height must be $h, so scale down to $h = k * LO
            $e->scale($h / $eh);
            my $ew_new = int(0.5 + $ew * $h / $eh);
            $p->_add_eps($e, $xbase + @$r * $w + 10, $ysize - ($ybase + $i*$h+$h));
        }
        else
        {
            $p->setfont("Times", 14);
            $p->setcolour("black");
            $p->text({align => "left", rotate => 0}, $xbase + @$r * $w + 10, $ysize - ($ybase + $i*$h+$h/2+4), $name);
        }
    }

    $p->setcolour(@$vcol);    
    $p->line($xbase + @$a_ta_H * $w, $ysize - $ybase, $xbase + @$a_ta_H * $w, $ysize - ($ybase + @$a_ta_M * $h));
    $p->line($xbase, $ysize - $ybase, $xbase + @$a_ta_H * $w, $ysize - $ybase);
    $p->line($xbase, $ysize - ($ybase + @$a_ta_M * $h), $xbase + @$a_ta_H * $w, $ysize - ($ybase + @$a_ta_M * $h));
    $p->line($xbase, $ysize - $ybase, $xbase , $ysize - ($ybase+ @$a_ta_M * $h));

    drawScale ($xbase/4, $ybase+$s_y, $min, $mid, $max, $sep, $res, $p, $xsize, $ysize, $scalefont, $s_h, $s_w, $uppertext, $lowertext, $lcol, $hcol, $mcol);

    my $outpdf = $fn ;
    my $outeps = $outpdf ;
    $outeps =~ s/pdf$/eps/ ;
    $p->output("$outeps") ;
    system("ps2pdf -dEPSCrop -dAutoRotatePages=/None $outeps $outpdf") ;
    system("rm -r ./Temp") ;

    unshift(@$a_ta_H, $dummy) ;
    unshift(@$a_ta_M, $a_ta_H) ;
}

 




sub printHeatmap
{
    my ($a_ta_M, $fn, $min, $mid, $max, $h, $w, $xplus, $yplus, $header_motif, $row_motif, $s_y, $s_h, $s_w, $sep, $res, $scalefont, $uppertext, $lowertext, $lcol, $hcol, $mcol) = @_ ;
    
    my $a_ta_H = shift(@$a_ta_M) ;
    my $dummy = shift(@$a_ta_H) ;

    my $xbase        = 100 ;
    my $ybase        = 150 ;
    my $ysize        = $ybase + $h * scalar(@$a_ta_M) + $yplus ;
    my $xsize        = $xbase + $w * scalar(@$a_ta_H) + $xplus ;
    
    my $p = new PostScript::Simple(xsize     => $xsize,
				   ysize     => $ysize,
				   colour    => 1,
				   eps       => 1,
				   units     => "pt") ;
    $p->setlinewidth(0.5) ;
    
    if ($header_motif eq "true" or $row_motif eq "true")
    {
	system("mkdir ./Temp") if !(-d "./Temp") ;
    }

    for (my $j=0; $j<@$a_ta_H; $j++) 
    {
	if ($header_motif eq "true")
	{
	    my $motif = $a_ta_H->[$j] ;
	    my $mo = Sets::myre2wm($motif);
	    open OUT, "> ./Temp/$motif.txt" or die "cannot open\n";
	    print OUT $mo;
	    close OUT;
	    system("../FIRE/SCRIPTS/weblogo/seqlogo -f ./Temp/$motif.txt -F EPS  -a -c -M -n -Y -w 5 -h 3 > ./Temp/$motif.eps");

	    my $e  = new PostScript::Simple::EPS(file => "./Temp/$motif.eps");
	    my $eh = $e->height;
	    my $ew = $e->width;

	    # height must be $h, so scale down to $h = k * LO
	    $e->scale($w / $eh);
	    $e->rotate(90);
	    my $ew_new = int(0.5 + $ew * $h / $eh);
	    $p->_add_eps($e, $xbase + ($j+1.1) * $w,  $ysize - ($ybase-3)); 
	}
	else
	{
	    $p->setcolour("black");
	    $p->setfont("Courrier", 6);
	    $p->text( { rotate => 90 }, $xbase + $j * $w + 3*$w/4, $ysize - ($ybase-3), $a_ta_H->[$j]);
	}
    }

    my @col = () ;
    for (my $i=0; $i<@$a_ta_M; $i++)
    {
	my $r  = $a_ta_M->[$i] ;
	my $name = shift @$r;
	print $name,"\n" ;
	for (my $j=0; $j<@$r; $j++)
	{
	    my $v = $r->[$j] ;

	    #defining the color
	    if (defined $mcol)
	    {
		if ($v>$mid)
		{
		    @col = Sets::interp_general( $v, $mcol, $hcol, $mid, $max);
		}
		else
		{
		    @col = Sets::interp_general( $v, $lcol, $mcol, $min, $mid);
		}
	    }
	    else
	    {
		@col = Sets::interp_general( $v, $lcol, $hcol, $min, $max);
	    }

	    $p->setcolour(@col);
	    $p->box({filled => 1},
		    $xbase + $j * $w,      $ysize - ($ybase + $i*$h) ,
		    $xbase + $j * $w + $w, $ysize - ($ybase + ($i*$h+$h)));
	}
			
	if ($row_motif eq "true")
        {
            my $motif = $name ;
            my $mo = Sets::myre2wm($motif);
            open OUT, "> ./Temp/$motif.txt" or die "cannot open\n";
            print OUT $mo;
            close OUT;
            system("../FIRE/SCRIPTS/weblogo/seqlogo -f ./Temp/$motif.txt -F EPS  -a -c -M -n -Y -w 5 -h 3 > ./Temp/$motif.eps");

            my $e  = new PostScript::Simple::EPS(file => "./Temp/$motif.eps");
            my $eh = $e->height;
            my $ew = $e->width;

            # height must be $h, so scale down to $h = k * LO
            $e->scale($h / $eh);
            my $ew_new = int(0.5 + $ew * $h / $eh);
            $p->_add_eps($e, $xbase + @$r * $w + 10, $ysize - ($ybase + $i*$h+$h));
        }
	else
	{
	    $p->setfont("Times", 14);
	    $p->setcolour("black");
	    $p->text({align => "left", rotate => 0}, $xbase + @$r * $w + 10, $ysize - ($ybase + $i*$h+$h/2+4), $name);
	}
    }

    drawScale ($xbase/4, $ybase+$s_y, $min, $mid, $max, $sep, $res, $p, $xsize, $ysize, $scalefont, $s_h, $s_w, $uppertext, $lowertext, $lcol, $hcol, $mcol) ;

    my $outpdf = $fn ;
    my $outeps = $outpdf ;
    $outeps =~ s/pdf$/eps/ ;
    $p->output("$outeps") ;
    system("ps2pdf -dEPSCrop -dAutoRotatePages=/None $outeps $outpdf") ;
    system("rm -r ./Temp") ;

    unshift(@$a_ta_H, $dummy) ;
    unshift(@$a_ta_M, $a_ta_H) ;
}

sub drawScale
{
    my ($x, $y, $min, $mid, $max, $sep, $res, $p, $xsize, $ysize, $scalefont, $h, $w, $uppertext, $lowertext, $lcol, $hcol, $mcol) = @_;
    
    $p->setcolour("black");
    $p->setfont("Courrier", $scalefont);
    $p->text({align => "center"}, $x+$w/2+0, $ysize - ($y - 3), $max);
    $p->text({align => "center"}, $x+$w/2+0, $ysize - ($y - 13), $uppertext);


    my $t = $max;
    for (my $i=0; $i<$res/2; $i++)
    {
	my @col = ();
	
	if (defined $mcol)
        {
            if ($t>$mid)
            {
                @col = Sets::interp_general( $t, $mcol, $hcol, $mid, $max);
            }
            else
            {
                @col = Sets::interp_general( $t, $lcol, $mcol, $min, $mid);
            }
        }
        else
        {
            @col = Sets::interp_general( $t, $lcol, $hcol, $min, $max);
	    $t -= ($max - $min) / $res;
        }

	$p->setcolour( @col );
	$p->box({filled => 1}, $x, $ysize - ($y + $i*$h) , $x+$w, $ysize - ($y  + $i*$h + $h));
	$t -= ($max - $min) / $res;
    }

    $p->setcolour("black");
    $p->setfont("Courrier", $scalefont);

    $p->text({align => "center"}, $x+$w/2+0, $ysize - ($y + ($res/2)*$h + $sep*3/4 ), $mid);

    if (defined $mcol)
    {    
	for (my $i=$res-1; $i>=$res/2; $i--)
	{
	    my @col = ();

	    @col = Sets::interp_general( $t, $mcol, $lcol, $min, $max);

	    $p->setcolour( @col );
	    $p->box({filled => 1}, $x, $ysize - ($y + $sep + $i*$h) , $x+$w, $ysize - ($y + $sep + $i*$h + $h));
	    $t -= ($max - $min) / $res;
	}     
	$p->setcolour( 'black' );           
	$p->text({align => "center"}, $x+$w/2+1, $ysize - ($y + $sep+$res*$h + 20 ), $lowertext);    
	$p->text({align => "center"}, $x+$w/2+0, $ysize - ($y + $sep+$res*$h + 10 ), $min);
    }
}


sub saveTable
{
    my ($a_ta_M, $fn) = @_ ;
    open O, "> $fn" or die "couldn't open $fn..." ;
    
    foreach my $r (@$a_ta_M)
    {
	print O join("\t", @$r) . "\n";
    }
    close O ;
}

sub clusterTable
{
    #$status = 'r' (rows) 'c' (columns) 'rc' (both)
    my ($a_ta_M, $status) = @_ ;
    my $a_new_M ;
    $a_new_M = &clusterTableRows($a_ta_M) if $status eq "r";
    $a_new_M = &clusterTableCols($a_ta_M) if $status eq "c";
    if ($status eq "rc" or $status eq "cr")
    {
	$a_new_M = &clusterTableRows($a_ta_M) ;
	$a_new_M = &clusterTableCols($a_new_M) ;
    }
    return $a_new_M ;
}

sub clusterTableRows
{
    my $ref_M = shift @_ ;
    my $a_ta_M = copyTable($ref_M) ;
    my $a_ta_H = shift(@$a_ta_M) ;
    my $dummy = shift(@$a_ta_H) ;
    my $a_ta_R ;
    foreach my $r (@$a_ta_M)
    {
	push(@$a_ta_R, $r->[0]) ;
	shift @$r ;
    }
    
    my %param = (data       => $a_ta_M,
		 mask       => '',
		 weight     => '',
		 transpose  => 0,
		 dist       => 'c',
		 method     => 'a') ;
    
    my ($result, $linkdist) = Algorithm::Cluster::treecluster(%param) ;
    die "The clustering failed..." if $result == undef ;

    my @a ;
    buildTree($result->[@$result-1], \@a, $result) ;
    
    my $new_M ;
    foreach my $i (@a)
    {
	push(@$new_M, $a_ta_M->[$i]) ;
	unshift(@{$new_M->[@$new_M-1]}, $a_ta_R->[$i]) ;
    }
    unshift (@$a_ta_H, $dummy) ;
    unshift (@$new_M, $a_ta_H) ;
 
    return $new_M ;
}


sub clusterTableCols
{
    my $ref_M = shift @_ ;
    my $a_ta_M = copyTable($ref_M) ;
    my $a_ta_H = shift(@$a_ta_M) ;
    my $dummy = shift(@$a_ta_H) ;
    my $a_ta_R ;
    foreach my $r (@$a_ta_M)
    {
        push(@$a_ta_R, $r->[0]) ;
        shift @$r ;
    }

    my %param = (data       => $a_ta_M,
                 mask       => '',
                 weight     => '',
                 transpose  => 1,
                 dist       => 'c',
                 method     => 'a') ;

    my ($result, $linkdist) = Algorithm::Cluster::treecluster(%param) ;
    die "The clustering failed..." if $result == undef ;
    
    my @a ;
    buildTree($result->[@$result-1], \@a, $result) ;

    my $new_M ;
    my $new_H ;
    my $c=0 ;
    foreach my $j (@a)
    {
	$new_H->[$c] = $a_ta_H->[$j] ;
	for (my $i=0 ; $i<@$a_ta_R ; $i++)
	{
	    $new_M->[$i]->[$c] = $a_ta_M->[$i]->[$j] ;
	}
	$c++ ;
    }
    for (my $i=0 ; $i<@$a_ta_R ; $i++)
    {
        unshift(@{$new_M->[$i]}, $a_ta_R->[$i]) ;
    }
    unshift (@$new_H, $dummy) ;
    unshift (@$new_M, $new_H) ;

    return $new_M ;
}

sub buildTree
{
    my ($n, $a, $t) = @_;
    push(@$a, $n->[0]) if ($n->[0]>=0) ;
    push(@$a, $n->[1]) if ($n->[1]>=0) ;
    buildTree($t->[-1*$n->[0]-1], $a, $t) if ($n->[0]<0) ;
    buildTree($t->[-1*$n->[1]-1], $a, $t) if ($n->[1]<0) ;
}

sub kmeansTable
{
    my $ref_M = shift @_ ;
    my $a_ta_M = copyTable($ref_M) ;
    my $a_ta_H = shift(@$a_ta_M) ;
    my $dummy = shift(@$a_ta_H) ;
    my $a_ta_R ;
    foreach my $r (@$a_ta_M)
    {
        push(@$a_ta_R, $r->[0]) ;
        shift @$r ;
    }
    my $mask ;
    for (my $i=0 ; $i<@$a_ta_R ; $i++)
    {
	for (my $j=0 ; $j<@$a_ta_H ; $j++)
	{
	    if ($a_ta_M->[$i]->[$j] eq "" or $a_ta_M->[$i]->[$j] eq "nan" or $a_ta_M->[$i]->[$j] eq "NaN" or $a_ta_M->[$i]->[$j] eq "NA")
	    {
		$mask->[$i]->[$j] = 0 ;
	    }
	    else
	    {
		$mask->[$i]->[$j] = 1 ;
	    }
	}
    }
    my $nbclusters = int(sqrt(scalar(@$a_ta_R))) ;
    print "Clustering ... $nbclusters clusters\n" ;

    my %param = (nclusters => $nbclusters,
		 data       => $a_ta_M,
                 mask       => '',
                 weight     => '',
                 transpose  => 0,
		 npass      => 10,
                 dist       => 'c',
		 method     => 'a',
		 initialid  => []) ;
    
    my ($clusters, $error, $found) = Algorithm::Cluster::kcluster(%param) ;
    
    my $new_M ;
    for (my $i=0 ; $i<@$a_ta_R ; $i++)
    {
	$new_M->[$i]->[0] = $a_ta_R->[$i] ;
	$new_M->[$i]->[1] = $clusters->[$i] ;
    }
    return $new_M ;
}

1 ;
