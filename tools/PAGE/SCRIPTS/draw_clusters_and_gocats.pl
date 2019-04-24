use lib "$ENV{PAGEDIR}/SCRIPTS";
use lib "$ENV{PAGEDIR}/SCRIPTS/PostScript-Simple-0.07/lib";

my $pagedir    = $ENV{PAGEDIR};
my $programdir = $pagedir."/PROGRAMS" ;
my $scriptdir  = $pagedir."/SCRIPTS" ;

use Getopt::Long;
use Fire;
use Table;
use Sets;
use PostScript::Simple;

use strict;

my $nbclusters = 0 ;
my $expfile     = undef;
my $data        = undef;
my $summaryfile = undef;
my $pvaluematrixfile = undef ;
my $ps2pdf      = 1;

GetOptions ('expfile=s'              => \$expfile,
	    'data=s'                 => \$data,
	    'pvaluematrixfile=s'     => \$pvaluematrixfile);


my $file = Sets::filename($expfile);
my $dir = "$expfile\_PAGE" ;
my $outclu = "$dir/clusters.txt";

my $ta = Table->new;


my %CLUSTERS = ();
open(EXP, "< $expfile") ;
<EXP> ;
while(<EXP>)
{
    s/\s+$// ;
    my($gene, $c) = split(/\t/, $_) ;
    push @{$CLUSTERS{$c}}, $gene;
}

system("perl -pi -e 's/\r//g' $data") ;
my %DATA = ();
$ta->loadFile($data);
my $a_ref_data = $ta->getArray();

my $A_REF_CONDS = shift @$a_ref_data;
shift @$A_REF_CONDS;
foreach my $r (@$a_ref_data) {
 my $n = shift @$r;
 $DATA{ $n } = $r;
}

my %PVALUES = undef;
if (defined($pvaluematrixfile)) {
 print "loading $pvaluematrixfile)\n" if (-e $pvaluematrixfile);
 my $h_ref_pv = $ta->getBidimensionalHash($pvaluematrixfile);
 %PVALUES = %$h_ref_pv;
}


#
# calculate the centroids
#
my $data_size = 0 ;
my %CENTROIDS = ();
foreach my $c (keys(%CLUSTERS)) {

 my $s        = [];
 my @cnt_rows = ();

 foreach my $r (@{$CLUSTERS{$c}}) {
   for (my $i=0; $i<@{ $DATA{$r} }; $i++) {
     if ($DATA{$r}->[$i] ne "") {
       $s->[$i] += $DATA{$r}->[$i];
       $cnt_rows[$i]++;
     }
   }
   # $s = Sets::addArrays($s, $DATA{$r});
   $data_size = @{$DATA{$r}} ;
 }

 my $n = @{$CLUSTERS{$c}};
 for (my $i=0; $i<@$s; $i++) {

   if ($cnt_rows[$i] == 0) {
     die "Problem, col $i has no values.\n";
   }

   $s->[$i] /= $cnt_rows[$i];
 }
 $CENTROIDS{$c} = $s;
}


my %SAMPLE_CEN ;
my %SAMPLE_NUM ;
my %SAMPLES ;
my $c ;
foreach $c (sort { $a <=> $b } (keys(%CENTROIDS)))
{
    for (my $i=0 ; $i<@$A_REF_CONDS ; $i++)
    {
	$SAMPLE_CEN{$c}{$A_REF_CONDS->[$i]}+=$CENTROIDS{$c}->[$i] ;
	$SAMPLE_NUM{$c}{$A_REF_CONDS->[$i]} ++ ;
	$SAMPLES{$A_REF_CONDS->[$i]} = 1 if (! (defined $SAMPLES{$A_REF_CONDS->[$i]})) ;
    }
}
foreach my $c (sort { $a <=> $b } (keys(%CENTROIDS)))
{
    foreach my $s (keys %SAMPLES)
    {
	$SAMPLE_CEN{$c}{$s}/=$SAMPLE_NUM{$c}{$s} ;
    }
}

#
#  print out centroids here
#
open OUTC, ">$outclu" or die "Cannot open $outclu.\n";
foreach my $c (sort { $a <=> $b } (keys(%CENTROIDS))) {
    print OUTC "$c\t" . join("\t", @{$CENTROIDS{$c}}) . "\n";
}

close OUTC; 

open OS, "> $dir/sample_centroids.txt" or die "Cannot open $outclu.\n";
print OS "Clusters" ;
foreach my $s (sort keys %SAMPLES)
{
    print OS "\t$s" ;
}
print OS "\n" ;
foreach my $c (sort { $a <=> $b } (keys(%CENTROIDS))) 
{
    print OS "$c" ;
    foreach my $s (sort keys %SAMPLES)
    {
	print OS "\t", $SAMPLE_CEN{$c}{$s} ;
    }
    print OS "\n" ;
}

close OS; 


