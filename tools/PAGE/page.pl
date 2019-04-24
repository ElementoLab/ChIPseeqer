#!/usr/bin/perl -w


my $pagedir ;
BEGIN{
    if ((!$ENV{PAGEDIR}) || ($ENV{PAGEDIR} eq '')) {
	$pagedir="./" ;
	print "The PAGEDIR environment variable is not set. It is set to default.\n";
    }
    else{
	$pagedir = $ENV{PAGEDIR};
    }
}

my $programdir = $pagedir."/PROGRAMS" ;
my $scriptdir  = $pagedir."/SCRIPTS" ;

use lib "$pagedir/SCRIPTS";


use Sets;
use PBS;
use Table;
use Getopt::Long;

use strict;


my $cmdline = "page.pl";
foreach my $r (@ARGV) {
  $cmdline .= " $r";
}

if (@ARGV == 0) {
  die "Usage: perl page.pl --expfile=FILE --datafile=FILE --goindexfile=FILE --gonamesfile=FILE --exptype=TXT --cattypes=F,C,P --catmaxcount=INT

Other options:
--ignore_non_signif=INT
--collabels=FILE

\n";
}


my $draw             = undef;
my $expfile          = undef ;
my $goindexfile      = undef ;
my $gonamesfile      = undef ;
my $species          = undef ;
my $exptype          = "discrete" ;
my $datafile         = undef ;
my $catmaxcount      = 300;
my $catmincount      = 0;
my $independence     = 1 ;
my $shuffle          = undef;
my $max_p            = 0.005 ;
my $randomize        = undef;
my $w                = undef;
my $minmax_lp        = undef;
my $minr             = 0;
my $nbclusters       = 5;
my $cattypes         = "P";  # undef means all (F,P,C)
my $verbose          = 0;
my $submit           = 0;
my $draw_sample_heatmap = "false" ;
my $draw_min         = -5;
my $draw_max         = 5;
my $ebins            = undef ; 
my $from_gui         = 0 ;
my $onecatfile       = undef;
my $suffix           = undef;
my $pathwaylist      = undef;
my $correl           = undef;
my $fontsize         = 20;
my $xsize            = undef;
my $outprofiles      = undef;
my $arrowcentertext  = undef;
my $arrowlefttext    = undef;
my $arrowrighttext   = undef;
my $font             = "Arial";
my $yscale           = 130;
my $dodraw           = 1;
my $dopage           = 1;
my $drawframes       = 0;
my $catstoignore     = undef;
my $collabels        = undef;
my $ybase            = undef;
my $onlysignif       = undef;
my $contheader       = undef;
my $colmap           = undef;
my $horizscale       = undef;
my $mirror           = undef;
my $pmax             = undef;
my $expdata          = undef;
my $ignore_non_signif= undef;
my $rotateheader     = undef;

GetOptions ('expfile=s'              => \$expfile,
	    'exptype=s'              => \$exptype,
	    'species=s'              => \$species,
	    'pathways=s'             => \$species,
	    "expdata=s"              => \$expdata,
	    'goindexfile=s'          => \$goindexfile,
	    "rotateheader=s"         => \$rotateheader,
	    'onecatfile=s'           => \$onecatfile,	    
	    'gonamesfile=s'          => \$gonamesfile,
	    'pathwaylist=s'          => \$pathwaylist,
	    'catstoignore=s'         => \$catstoignore,
	    'catmaxcount=s'          => \$catmaxcount,
	    'catmincount=s'          => \$catmincount,
	    'cattypes=s'             => \$cattypes,
	    'collabels=s'            => \$collabels,
	    'shuffle=s'              => \$shuffle,
	    'correl=s'               => \$correl,
	    'datafile=s'             => \$datafile,
	    'max_p=s'                => \$max_p,
	    'minmax_lp=s'            => \$minmax_lp,
	    'xsize=s'                => \$xsize,
	    'minr=s'                 => \$minr,
	    'font=s'                 => \$font,
	    'nbclusters=s'           => \$nbclusters,
	    'arrowcentertext=s'      => \$arrowcentertext,
	    'arrowlefttext=s'        => \$arrowlefttext,
	    'arrowrighttext=s'       => \$arrowrighttext,	    
	    'yscale=s'               => \$yscale,
	    'ybase=s'                => \$ybase,
	    'dopage=s'               => \$dopage,
	    'dodraw=s'               => \$dodraw,
	    'verbose=s'              => \$verbose,
	    'independence=s'         => \$independence,
	    'submit=s'               => \$submit,
	    'fontsize=s'             => \$fontsize,
	    'draw_sample_heatmap'    => \$draw_sample_heatmap,
	    'randomize=s'            => \$randomize,
	    'draw_min=s'             => \$draw_min,
	    'draw_max=s'             => \$draw_max,
	    'draw=s'                 => \$draw,
	    'ebins=s'                => \$ebins,
	    'suffix=s'               => \$suffix,
	    'outprofiles=s'          => \$outprofiles,
	    'from_gui=i'             => \$from_gui,
	    'onlysignif=s'           => \$onlysignif,
	    'contheader=s'           => \$contheader,
	    'colmap=s'               => \$colmap,
	    'horizscale=s'           => \$horizscale,
	    'mirror=s'               => \$mirror,
	    "pmax=s"                 => \$pmax,
	    "ignore_non_signif=s"    => \$ignore_non_signif,
	    "w=s"                    => \$w) ;

if (defined($suffix)) {
  
  if (($expfile =~ /\*/) or ($submit == 1)) {
    die "--suffix option not supported with multiple expfile, or grid submission (for now).\n";
  } else {

    if ($suffix =~ /^pathway/) {
      $suffix = $species;

      if ($cattypes ne "P") {
	$suffix .= ".$cattypes";
      }

    }

    my $newexpfile = "$expfile.$suffix";
    system("cp $expfile $newexpfile");
    $expfile = $newexpfile;
    
  }

}  

#
# Make a batch in case there are multiple input files
#
if (($expfile =~ /\*/) or ($submit == 1))
{
    my $files = Sets::getFiles($expfile) ;

    my $walltime = "20:00:00";
    my $platform = undef;
    
    foreach my $file(@$files)
    {
	my $f = substr($file, rindex($file, "/"), length($file)-rindex($file, "/")) ;
	mkdir "$file\_PAGE/" if ! (-d "$file\_PAGE/");
	my $expfile_nodups_page  = "$file\_PAGE/$f";
	
	my $pwd  = `pwd`; $pwd =~ s/\n//;
	my $time = Sets::getNiceDateTime(1);
	
	my $pbs = PBS->new;
	$pbs->setPlatform($platform) if (defined($platform));
	$pbs->setWallTime($walltime);
	$pbs->addCmd("cd $pwd");
	
	$pbs->setScriptName("$expfile_nodups_page.script");
	
	$pbs->addCmd("date") ;
	$pbs->addCmd("export PAGEDIR=$pwd") ;
	$pbs->addCmd("echo \"Running PAGE\"") ;
	
	my $cmd = "perl $pagedir/page.pl --expfile=$file --species=$species --goindexfile=$goindexfile --gonamesfile=$gonamesfile --exptype=$exptype --catmaxcount=$catmaxcount --cattypes=$cattypes --independence=$independence --datafile=$datafile --shuffle=$shuffle --max_p=$max_p --ebins=$ebins --draw_min=$draw_min --draw_max=$draw_max" ;
	$pbs->addCmd($cmd) ;
	
	my $page_jobid ;
	if ($submit==0)
	{
	    $pbs->execute ;
	}
	elsif ($submit==1)
	{
	    $page_jobid = $pbs->submit ;
	    print "Submitted job $page_jobid.\n";
	}
    }
    exit (1) ;
}

#
# Making workspace directory
#
if (! -d "$expfile\_PAGE") 
{
    mkdir("$expfile\_PAGE") or die "couldn't make the directory: $?";
}


#
# save command line to _FIRE
# 
open OUTC, ">$expfile\_PAGE/cmdline.txt" or print "Cannot open $expfile\_PAGE/cmdline.txt";
print OUTC "$cmdline\n";
close OUTC;



#
# Changing unconventional inputs
#
if (defined $species and $species ne "" and !defined($goindexfile)){
	if ($species =~ /kegg/i){
		# print STDERR "KEGG\n";
		$species= "kegg";
	}
  $goindexfile = "$pagedir/PAGE_DATA/ANNOTATIONS/$species/$species\_index.txt" ;
  $gonamesfile = "$pagedir/PAGE_DATA/ANNOTATIONS/$species/$species\_names.txt" ;
} elsif (defined($onecatfile)) {
  
  if (! -e $onecatfile) {
    die "Cannot open --onecatfile=$onecatfile\n";
  }

  $goindexfile = $onecatfile;
  my $f = "$onecatfile.names";
  #if (! -e $f) {
    
    my %ONECATNAMES = ();
    open IN, $onecatfile;
    my $l = <IN>;
    while (my $l = <IN>) {
      chomp $l;
      my @a = split /\t/, $l, -1;
      $ONECATNAMES{$a[1]} = 1;
    }
    close IN;



    # create a fake name file
    open OUT, ">$f" or die "Cannot open $f\n";
    foreach my $v (sort(keys(%ONECATNAMES))) {
      print OUT "$v\t$v\tP\n";
    }
    close OUT;
    

  #}

  
  $gonamesfile = $f;

  # remove limit for cat max count
  $catmaxcount = 10000000;
  
}



if ($exptype eq "discrete") {
  $exptype = 1 ;
} elsif ($exptype eq "continuous"){
  $exptype = 0 ;
}

if ($catmaxcount eq "all") {
  print INF "Retaining all categories.\n" ;
  $catmaxcount = -1 if ($catmaxcount eq "all") ;
}

my %PARAMS = (expfile          => $expfile,
	      goindexfile      => $goindexfile,
	      gonamesfile      => $gonamesfile,
	      catmaxcount      => $catmaxcount,
	      catmincount      => $catmincount,
	      pathwaylist      => $pathwaylist,
	      outprofiles      => $outprofiles,
	      exptype          => $exptype,
	      cattypes         => $cattypes,
	      max_p            => $max_p,
	      correl           => $correl,
	      verbose          => $verbose,
	      minr             => $minr,
	      independence     => $independence,
	      ebins            => $ebins);

#
# Write a log of the input data
#
open(INF, "> $expfile\_PAGE/info.txt") ;
foreach my $k (sort keys %PARAMS){

  print INF $k, "\t", $PARAMS{$k}, "\n" if (defined($PARAMS{$k}));
}

my $todo = &getPAGECommand(\%PARAMS);

if ($verbose == 1) {
  print "$todo\n";
}

if (($dopage == 1) && (!defined($randomize))) {
  print "$todo\n" if ($verbose == 1);
  system("$todo") == 0 or die "command failed: $todo";
}
  
#
# Randomizing the data
#
if (defined($randomize) && ($randomize > 0)) {
  
  my $rand_dir = "$expfile\_PAGE/RANDOM";
  if (! -e $rand_dir) {
    mkdir $rand_dir;
  }

  my @rands = ();
  for(my $i=0 ; $i<$randomize ; $i++) {
    
    my $rrun = $i+1;
    print "Randomize run $rrun. ";

    

    # create name
    my $rand_expfile = "$expfile\_PAGE/RANDOM/$i.txt";

    # PAGE dir
    my $rand_page_dir = "$rand_expfile\_PAGE";
    mkdir $rand_page_dir if (! -e $rand_page_dir);

    # create file
    system("perl $ENV{PAGEDIR}/SCRIPTS/shuffle_column.pl $expfile > $rand_expfile") == 0 or die "system failed: $?";

    # update parameter
    $PARAMS{expfile}       = $rand_expfile;
    $PARAMS{stdoutlogfile} = "$rand_page_dir/log.txt";

    # get cmd 
    my $todo = &getPAGECommand(\%PARAMS);
    if ($verbose == 1) {
      print "$todo\n";
    }
    system($todo) == 0 or die "system failed: $?";

    # parse output
    open(IN, "< $rand_page_dir/log.txt") or die "couldn't open file: $?" ;
    my $fp = undef;
    while(<IN>) {
      chomp ;
      if (/Number of categories that passed the tests/) {
	my @A = split(/\s+/, $_) ;
	#print("Number of false positives = ", $A[-1],"\n") ;
	$fp = $A[-1];
	push(@rands, $A[-1]) ;
      }
    }
    close(IN) ;

    print "Random run yielded #FP = $fp\n";
    
  }

  my $average = Sets::average(\@rands);
  print ("Results from $randomize random shufflings indicate a mean of $average false positives");
  if ($randomize > 1) {    
    my $std     = Sets::stddev (\@rands);
    print " with standard deviation of $std.\n";
  } else {
    print ".\n";
  }

  #foreach(@rands) {
  #  print INF "Number of false positives = ".$_."\n" ;
  #}
  #print INF "Results from $shuffle random shufflings indicate a mean of $average false positivies with standard deviation of $std\n" ;
}

close INF;


if (($dodraw == 1) && !defined($randomize))  {
  #
  # drawing the result
  #
  my $pvaluematrixfile = "$expfile\_PAGE/pvmatrix.txt";


  $todo = "perl $scriptdir/mi_go_draw_matrix.pl  --pvaluematrixfile=$pvaluematrixfile --expfile=$expfile --draw_sample_heatmap=$draw_sample_heatmap --min=$draw_min --max=$draw_max --cluster=$nbclusters" ; 
  
  if (defined($onecatfile)) {
    $todo .= " --onecat=1 ";
  }
  
  if (defined($ignore_non_signif)) {
    $todo .= " --ignore_non_signif=$ignore_non_signif ";
  }
  
  if (defined($minmax_lp)) {
    $todo .= " --minmax_lp=$minmax_lp ";
  }
  
  if ($exptype == 0) {
    $todo .= " --quantized=0 ";
  }
  
  if (defined($expdata)) {
    $todo .= " --expdata=$expdata ";
  }
  
  if (defined($fontsize)) {
    $todo .= " --fontsize=$fontsize";
  }
  if (defined($xsize)) {
    $todo .= " --xsize=$xsize";
  }
  
  if (defined($arrowcentertext)) {
    $todo .= " --arrowcentertext=\"$arrowcentertext\" ";
  }
  
  if (defined($arrowlefttext)) {
    $todo .= " --arrowlefttext=\"$arrowlefttext\" ";
  }
  
  if (defined($arrowrighttext)) {
    $todo .= " --arrowrighttext=\"$arrowrighttext\" ";
  }
 
  if (defined($rotateheader)) {
    $todo .= " --rotateheader=$rotateheader ";
  }
  
  if (defined($font)) {
    $todo .= " --font=$font ";
  }
  
  if (defined($yscale)) {
    $todo .= " --yscale=$yscale ";
  }

   if (defined($ybase)) {
    $todo .= " --ybase=$ybase ";
  }

  if (defined($catstoignore)) {
    $todo .= " --catstoignore=$catstoignore ";
  }
  
  if (defined($collabels)) {
    $todo .= " --collabels=$collabels ";
  }

  if (defined($onlysignif)) {
    $todo .= " --onlysignif=$onlysignif ";
  }

  if (defined($contheader)) {
    $todo .= " --contheader=$contheader ";
  }

  if (defined($colmap)) {
    $todo .= " --colmap=$colmap ";
  }
  
  if (defined($horizscale)) {
    $todo .= " --horizscale=$horizscale ";
  }
  
  if (defined($mirror)) {
    $todo .= " --mirror=$mirror ";
  }

  if (defined($pmax)) {
    $todo .= " --pmax=$pmax ";
  }

  if (defined($w)) {
    $todo .= " --w=$w ";
  }


  if ($verbose == 1) {
    print "$todo\n";
  }
  
  
  system("$todo") if ($from_gui==0) ;

  $todo = "perl $scriptdir/list_killed_cats.pl --pvmatrixfile=$pvaluematrixfile"; 
}




if (defined($draw)) {
  my $file   = Sets::filename($expfile);
  my $outpdf = "$expfile\_PAGE/$file.summary.pdf";
  system("$draw $outpdf");
}


sub getPAGECommand {
  my ($p) = @_;
  
  my $pvaluematrixfile = "$p->{expfile}\_PAGE/pvmatrix.txt";
  my $countmatrixfile  = "$p->{expfile}\_PAGE/countmatrix.txt";
  my $outquantized     = "$p->{expfile}\_PAGE/outquantized.txt";

  my $todo = "$programdir/page -expfile $p->{expfile} -goindexfile $p->{goindexfile} -gonamesfile $p->{gonamesfile} -catmaxcount $p->{catmaxcount} -catmincount $p->{catmincount} -logfile $pvaluematrixfile.log  -quantized $p->{exptype} -pvaluematrixfile $pvaluematrixfile -countmatrixfile $countmatrixfile -max_p $p->{max_p} ";
  
  if ($p->{exptype} == 0) {
    $todo .= " -outquantized $outquantized ";
  }

  my @a_cats = split /\,/, $p->{cattypes};
  if (Sets::in_array('F', @a_cats)) {
    $todo .= " -F 1 ";
  } else {
    $todo .= " -F 0 ";
  }

  if (Sets::in_array('P', @a_cats)) {
    $todo .= " -P 1 ";
  } else {
    $todo .= " -P 0 ";
  }

  if (Sets::in_array('C', @a_cats)) {
    $todo .= " -C 1 ";
  } else {
    $todo .= " -C 0 ";
  }

  if (defined($p->{minr})) {
    $todo .= " -minr $p->{minr} ";
  }
  
  if (defined($p->{correl})) {
    $todo .= " -correl $p->{correl} ";
  }

  if (defined($p->{verbose})) {
    $todo .= " -verbose $p->{verbose} ";
  }
  
  if (defined($p->{pathwaylist})) {
    $todo .= " -pathwaylist $p->{pathwaylist} ";
  }

  if (defined($p->{outprofiles})) {
    my $profile_filename = "$p->{expfile}\_PAGE/signif.profiles.txt";
    $todo .= " -outprofiles $profile_filename ";
  }

  if (defined($p->{stdoutlogfile})) {
    $todo .= " > $p->{stdoutlogfile} ";
  }

  if (defined($p->{independence})) {
    $todo .= " -independence $p->{independence} ";
  }

  if (defined($p->{ebins})) {
    $todo .= " -ebins $p->{ebins} "; 
  }
  
  return $todo;
}
