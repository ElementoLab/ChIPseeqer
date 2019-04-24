#
# input : CG number
# what  : fetch all 3'UTR, align using dialign
#


use lib qw(/home/olly/PERL_MODULES);

use Table;
use Sets;
use Sequence;
use Fasta;
use strict;
use Repeats;
#use Data::Dumper;

my $fa = Fasta->new;

my $se = Sequence->new;


my $ta = Table->new;

#
#  load the orthologs
#
$ta->loadFile("DATA/ortholog_table.txt.6only");
my $a_ref_orth_info = $ta->getArray();
my @a_orthologs = ();
my $nb = 1;
foreach my $r (@$a_ref_orth_info) {

    next if (($r->[0] eq "agam") || ($r->[0] eq "amel"));
    
    $ta->loadFile("DATA/$r->[2]");
    
    my $h_ref_o = $ta->getIndex(0);

    my %h_tmp = (ORTHOLOGS => $h_ref_o, FILE => $r->[1], NAME => $r->[0], NUMBER => $nb);
    
    push @a_orthologs, \%h_tmp;

    $nb ++;
}




#
#  load the exons
#

$ta->loadFile("DATA/ensembl_droso_exons_CDS_reordered.tsv");
my $a_ref = $ta->getArray();

my %CHR            = ();
my %EXONS          = ();
my %BOUNDARIES     = ();
my %BOUNDARIES_CDS = ();
my %STRAND         = ();

foreach my $r (@$a_ref) {
    my @a_exon = ($r->[2], $r->[3]);
    push @{ $EXONS{ $r->[1] } }, \@a_exon;
    $CHR{ $r->[1] }             = $r->[0]; 
    $BOUNDARIES{ $r->[1] }->[0] = (defined($BOUNDARIES{ $r->[1] }->[0])?Sets::min($BOUNDARIES{ $r->[1] }->[0], $r->[2]):$r->[2]);
    $BOUNDARIES{ $r->[1] }->[1] = (defined($BOUNDARIES{ $r->[1] }->[1])?Sets::max($BOUNDARIES{ $r->[1] }->[1], $r->[3]):$r->[3]);
    $STRAND{ $r->[1] }          = $r->[4];

    # no need to do that if this exon is not part of the CDS
    next if (($r->[5] eq "") && ($r->[6] eq ""));

    #print Dumper($r);
    
    $BOUNDARIES_CDS{ $r->[1] }->[0] = (defined($BOUNDARIES_CDS{ $r->[1] }->[0])?Sets::min($BOUNDARIES_CDS{ $r->[1] }->[0], $r->[5]):$r->[5]);
    $BOUNDARIES_CDS{ $r->[1] }->[1] = (defined($BOUNDARIES_CDS{ $r->[1] }->[1])?Sets::max($BOUNDARIES_CDS{ $r->[1] }->[1], $r->[6]):$r->[6]);
}



# 
#  get the 3'UTR
#
my $start = undef;
my $end   = undef;
my $len   = undef;
my $seq   = undef;

if ($STRAND{ $ARGV[0] } == 1) {
    $start = $BOUNDARIES_CDS{ $ARGV[0] }->[1];
    $end   = $BOUNDARIES    { $ARGV[0] }->[1]; 
    $se->setBlastDB("DATA/dmel.fasta");
    $seq = $se->getSequenceFromBlastDB($CHR{ $ARGV[0] }, $start, $end);
    $len = length($seq);
} else {
    $start = $BOUNDARIES    { $ARGV[0] }->[0];
    $end   = $BOUNDARIES_CDS{ $ARGV[0] }->[0]; 
    $se->setBlastDB("DATA/dmel.fasta");
    $seq = $se->getSequenceFromBlastDB($CHR{ $ARGV[0] }, $start, $end);
    $seq = Sets::getComplement($seq);
    $len = length($seq);
}

exit if (length($seq) < 10);

open OUT, ">$ARGV[0].seq";

print OUT ">dmel\n$seq\n";

#
#  get the orthologous 3'UTRs, but larger ..
#
foreach my $o (@a_orthologs) {
    next if (!defined($o->{ORTHOLOGS}->{ $ARGV[0] }));
    
    die "oh oh\n" if ($o->{ORTHOLOGS}->{ $ARGV[0] }->[ 3 ] < $o->{ORTHOLOGS}->{ $ARGV[0] }->[ 2 ]);

    $se->setBlastDB("DATA/$o->{FILE}");
    my $st = undef;
    my $en = undef;

    my $ch = $o->{ORTHOLOGS}->{ $ARGV[0] }->[ 1 ];
    my $seq_orth = undef;
    if ($o->{ORTHOLOGS}->{ $ARGV[0] }->[ 4 ] == 1) {
	$st = $o->{ORTHOLOGS}->{ $ARGV[0] }->[ 3 ] - 200; 
	$st = 0 if ($st < 0);
	$en = $o->{ORTHOLOGS}->{ $ARGV[0] }->[ 3 ] + 200 + $len;	
	$seq_orth = $se->getSequenceFromBlastDB($ch, $st, $en); 
    } else {
	$st = $o->{ORTHOLOGS}->{ $ARGV[0] }->[ 2 ] - 200 - $len;
	$st = 0 if ($st < 0);
	$en = $o->{ORTHOLOGS}->{ $ARGV[0] }->[ 2 ] + 200;
	$seq_orth = $se->getSequenceFromBlastDB($ch, $st, $en); 
	$seq_orth = Sets::getComplement($seq_orth);
    }
    if (defined($seq_orth)) {
	print OUT ">$o->{NAME}\n$seq_orth\n";
    }
    
} 
   

close OUT;


my $todo = "~/RPMS/dialign2_dir/dialign2-2 -n $ARGV[0].seq";

system($todo);


