BEGIN{
    if ((!$ENV{PAGEDIR}) || ($ENV{PAGEDIR} eq '')) {
        print "The PAGEDIR environment variable is not set. Please set it using export or setenv (see PAGE tutorial online).\n";
        exit;
    }
}

my  $pagedir = $ENV{PAGEDIR};
my $programdir = $pagedir."/PROGRAMS" ;
my $scriptdir  = $pagedir."/SCRIPTS" ;

use lib "$ENV{PAGEDIR}/SCRIPTS";

use strict;
use Sets;
use Subs;
use PBS;
use Table;
use Getopt::Long;
use Data::Dumper;
use FileHandle ;

if (@ARGV == 0) 
{
    die "Usage: perl prmg.pl --expfile=FILE --firefile=FILE --pagefile=FILE --goindexfile=FILE --max_p=F (--tu_enabled=0/1 --tu_file=FILE --printmode=0/1) --dorna=0/1\n" ;
}

my $expfile     = undef ;
my $firefile    = undef ;
my $pagefile    = undef ;
my $goindexfile = undef ;
my $species     = undef ;
my $gonamesfile = undef ;
my $max_p       = 0.05 ;
my $submit      = 0 ;
my $printmode   = 0 ;
my $tu_enabled  = 0 ;
my $tu_file     = undef ;
my $dorna       = 1 ;
GetOptions ('expfile=s'     => \$expfile,
	    'firefile=s'    => \$firefile,
	    'pagefile=s'    => \$pagefile,
	    'goindexfile=s' => \$goindexfile,
	    'species=s'     => \$species,
	    'max_p=s'       => \$max_p,
	    'submit=s'      => \$submit,
	    'tu_enabled=s'  => \$tu_enabled,
	    'tu_file=s'     => \$tu_file,
	    'printmode=s'   => \$printmode,
	    'dorna=s'       => \$dorna) ;


if (($expfile =~ /\*/) or ($submit == 1))
{
    my $files = Sets::getFiles($expfile) ;

    my $walltime = "20:00:00";
    my $platform = undef;
    
    foreach my $file(@$files)
    {
	my $f = Sets::filename($file) ;
	mkdir "$file\_PAGE/" if ! (-d "$file\_PAGE/");
	my $expfile_nodups_prmg  = "$file\_PAGE/$f";
	
	my $pwd  = `pwd`; $pwd =~ s/\n//;
	my $time = Sets::getNiceDateTime(1);
	
	my $pbs = PBS->new;
	$pbs->setPlatform($platform) if (defined($platform));
	$pbs->setWallTime($walltime);
	$pbs->addCmd("cd $pwd");
	
	$pbs->setScriptName("$expfile_nodups_prmg.script");
	
	$pbs->addCmd("date") ;
	$pbs->addCmd("echo \"Running PRMG\"") ;
	
	my $firefile=$file."_FIRE/DNA_RNA/$f.summary" ;
	my $pagefile=$file."_PAGE/pvmatrix.txt" ;
	
	my $cmd = "prmg.pl --expfile=$file --firefile=$firefile --pagefile=$pagefile --goindexfile=$goindexfile --max_p=$max_p --dorna=$dorna" ;
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

if ($printmode == 1)
{
    my $dir = Sets::dirname($pagefile) ;

    my $ta=Table->new ;
    $ta->loadFile("$dir/motif_cat.cdt") ;
    my $a_ta_M = $ta->getArray() ;
    
    #$my ($a_ta_M, $fn, $vline, $hline, $vcol, $min, $mid, $max, $h, $w, $xplus, $yplus, $header_motif, $row_motif, $s_y, $s_h, $s_w, $sep, $res, $scalefont,$uppertext, $lowertext,$lcol, $hcol, $mcol) = @_ ;
    Subs::printHeatmapBars($a_ta_M, "$dir/motif_cat.pdf", "true","true",[0,0,0],-3, 0,3, 35, 30, 560, 50, "true", "false", 0, 2, 20, 8, 30, 12, "Pos", "Neg", [0,0,255], [255,0,0], [255,255,255]);

    die ;
}

if (defined $species and $species ne "" and $goindexfile == undef){
    $goindexfile = "$pagedir/PAGE_DATA/ANNOTATIONS/$species/$species\_index.txt" ;
}
if (! defined ($firefile) and $dorna==1){
    my $fn = Sets::filename($expfile) ;
    $firefile = "$expfile\_FIRE/DNA_RNA/$fn.summary" ;
}
if (! defined ($firefile) and $dorna==0){
    my $fn = Sets::filename($expfile) ;
    $firefile = "$expfile\_FIRE/DNA/$fn.summary" ;
}
if (! defined ($pagefile)){
    $pagefile = "$expfile\_PAGE/pvmatrix.txt" ;
}

my %TU ;
if ($tu_enabled){
    open I, "< $tu_file" or die;
    while(<I>){
	s/\s+$// ;
	my ($tu, @a) = split(/\t/, $_) ;
	$TU{$tu} = \@a ;
    }
}

#
# Loading expfile
#
print "Loading expfile ..." ;
my $ta = Table->new ;
$ta->loadFile($expfile) ;
my $a_expfile = $ta->getArray() ;
shift @$a_expfile ;
print "Done\n" ;

#
# Loading pvmatrix file
#
print "Loading pagefile ..." ;
$ta->loadFile($pagefile) ;
my $a_pv_page = $ta->getArray() ;
shift @$a_pv_page ;
my $a_catlist ;
my $c=0 ;
foreach my $r (@$a_pv_page)
{
    $r->[0] =~ s/,//;
    my (@f) = split(/\s/, $r->[0]) ;
    my $go = shift(@f) ;
    my $cat = join(" ", @f) ;
    $a_catlist->[$c]->[0] = $go ;
    $a_catlist->[$c]->[1] = $cat ;
    if ($cat eq ""){
	$a_catlist->[$c]->[1] = $go ;
    }
    $a_catlist->[$c]->[2] = "P" ;
    $c++ ;
}
Subs::saveTable($a_catlist, "$pagedir/TEMP/go_names.in") ;
print "Done\n" ;

#
# Loading motifs
#
$ta->loadFile("$firefile") ;
my @motif = @{$ta->getColumn(0)} ;

my $firebinfile = substr($firefile, 0, rindex($firefile, "DNA_RNA")) ;
my $filename = Sets::filename($firefile) ;
$filename =~ s/\.summary/\.profiles/ ;
$ta->loadFile("$firebinfile/DNA/$filename");
my $a_dna_pro = $ta->getArray() ;
my %motif_gene ;
my %gene_motif ;
print "Loading fire profile for DNA..." ;
foreach my $r (@$a_dna_pro)
{
    push(@{$motif_gene{$r->[0]}}, $r->[1]) ;
    push(@{$gene_motif{$r->[1]}}, $r->[0]) ;
    #push(@motif, $r->[0]) if ! (grep {$_ eq $r->[0]} (@motif)) ;
}
print "Done\n" ;
$ta->loadFile("$firebinfile/RNA/$filename");
my $a_rna_pro = $ta->getArray() ;
if ($dorna==1){
    print "Loading fire profile for RNA..." ;
    foreach my $r (@$a_rna_pro)
    {
	$r->[0] =~ s/T/U/g ;
	push(@{$motif_gene{$r->[0]}}, $r->[1]) ;
	push(@{$gene_motif{$r->[1]}}, $r->[0]) ;
	#push(@motif, $r->[0]) if ! (grep {$_ eq $r->[0]} (@motif)) ;
    }
    print "Done\n" ;
}

print "Loading motif names for DNA..." ;
my %motif_names ;
$filename =~ s/profiles$/motifnames/ ;
$ta->loadFile("$firebinfile/DNA/$filename");
my $a_motif_name = $ta->getArray() ;
foreach my $r (@$a_motif_name)
{
    my $motif = $r->[0] ;
    my $name = $r->[1] ;
    $name =~ s/^J\_// ;
    $name =~ s/^M\S(\d+)\_// ;
    $name =~ s/\.txt$// ;
    $motif_names{$motif} = $name ;
}
print "Done\n" ;
if ($dorna==1){
    print "Loading motif names for RNA..." ;
    $ta->loadFile("$firebinfile/RNA/$filename");
    $a_motif_name = $ta->getArray() ;
    foreach my $r (@$a_motif_name)
    {
	my $motif = $r->[0] ;
	my $name = $r->[1] ;
	$name =~ s/^J\_// ;
	$name =~ s/^M\S(\d+)\_// ;
	$name =~ s/\.txt$// ;
	$motif_names{$motif} = $name ;
    }
    print "Done\n" ;
}

my $GOINDEX = FileHandle->new("< $goindexfile") ;
my %go_gene ;
my %gene_go ;
while(<$GOINDEX>)
{
    s/\s+$// ;
    my ($gene, @a) = split(/\t/, $_) ;
    
    foreach my $i (0..$#a)
    {
	push(@{$go_gene{$a[$i]}}, $gene) ;
	push(@{$gene_go{$gene}}, $a[$i]) ;
    }
}

print "Running...\n" ;
my %h_pv ;
my $count=0 ;
foreach my $motif (@motif)
{
    print $motif , " ($count/", $#motif+1, ")";
    my $a_mo_pro ;
    foreach my $r (@$a_expfile)
    {
	my $g = $r->[0] ;
	if (grep {$motif eq $_} (@{$gene_motif{$g}}))
	{
	    if ($tu_enabled){
		next if (! defined $TU{$g}) ;
		my @genes = @{$TU{$g}} ;
		foreach my $b (@genes){
		    push(@$a_mo_pro, [$b, 1]) ;
		}
	    }
	    else{
		push(@$a_mo_pro, [$g, 1]) ;
	    }
	}
	else
	{
	    if ($tu_enabled){
		next if (! defined $TU{$g}) ;
                my @genes = @{$TU{$g}} ;
                foreach my $b (@genes){
                    push(@$a_mo_pro, [$b, 0]) ;
                }
            }
            else{
                push(@$a_mo_pro, [$g, 0]) ;
            }
	}
    }
    Subs::saveTable($a_mo_pro, "$pagedir/TEMP/motif_profile.in") ;

    my $cmd = "$pagedir/PROGRAMS/mi_go_motif_calculator -expfile $pagedir/TEMP/motif_profile.in -gonamesfile $pagedir/TEMP/go_names.in -goindexfile $goindexfile -catmaxcount -1 -P 1 -quantized 1 -independence 0 -max_p $max_p -outfile $pagedir/TEMP/temp" ;
    system ($cmd) ;

    open MI, "< $pagedir/TEMP/temp" ;
    while(<MI>)
    {
	chomp ;
	my ($go, $pvalue) = split(/\t/, $_) ;
	$h_pv{$go}{$motif} = $pvalue ;
    }
    print "\n" ;
    $count++ ;
}

open M, "> $pagedir/TEMP/temp" ;
foreach my $motif (@motif)
{
    print M "\t$motif" ;
    if (defined $motif_names{$motif} and $motif_names{$motif} ne "-")
    {
	print M "/", $motif_names{$motif} ;
    }
}
print M "\n" ;
foreach my $r (@$a_catlist)
{
    print M $r->[1]," ", $r->[0] ;
    print $r->[1]," ", $r->[0], "\n" ;

    foreach my $motif (@motif)
    {
	print M "\t", $h_pv{$r->[1]}{$motif} if (defined $h_pv{$r->[1]}{$motif}) ;
	print M "\t0" if !(defined $h_pv{$r->[1]}{$motif}) ;
    }
    print M "\n" ;
}
close M ;

$ta->loadFile("$pagedir/TEMP/temp") ;
my $a_ta_M = $ta->getArray() ;
$a_ta_M = Subs::clusterTable($a_ta_M, "rc") ;

my $dir = Sets::dirname($pagefile) ;
Subs::saveTable($a_ta_M, "$dir/motif_cat.cdt") ;

my $ta=Table->new ;
$ta->loadFile("$dir/motif_cat.cdt") ;
my $a_ta_M = $ta->getArray() ;

#$my ($a_ta_M, $fn, $vline, $hline, $vcol, $min, $mid, $max, $h, $w, $xplus, $yplus, $header_motif, $row_motif, $s_y, $s_h, $s_w, $sep, $res, $scalefont,$uppertext, $lowertext,$lcol, $hcol, $mcol) = @_ ;
Subs::printHeatmapBars($a_ta_M, "$dir/motif_cat.pdf", "true","true",[0,0,0],-3, 0,3, 35, 30, 560, 50, "true", "false", 0, 2, 20, 8, 30, 12, "Pos", "Neg", [0,0,255], [255,0,0], [255,255,255]);
