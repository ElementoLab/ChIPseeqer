#
# input : CG number
# what  : fetch all 3'UTR, align using dialign
#


BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";



use Table;
use Sets;
use Sequence;
use Fasta;
use strict;
use DataFiles;

my $verbose = 1;

my $df = DataFiles->new;

my $fa = Fasta->new;

my $se = Sequence->new;
$se->setVerbose($verbose);

my $ta = Table->new;

if (scalar(@ARGV) == 0) {
  die "usage : gene_name [UD] len ortholog_summary\n";
}

my $len    = $ARGV[2];
my $offset = 200;

#
#  load the orthologs
#
$ta->loadFile($ARGV[3]);
my $a_ref_orth_info = $ta->getArray();

# get the reference species
my $a_ref_ref = shift @$a_ref_orth_info;


my @a_orthologs = ();
my $nb = 1;
foreach my $r (@$a_ref_orth_info) {    
    $ta->loadFile("DATA/$r->[2]");    
    my $h_ref_o = $ta->getIndex(0);
    my %h_tmp = (ORTHOLOGS => $h_ref_o, FILE => $r->[1], NAME => $r->[0], NUMBER => $nb);    
    push @a_orthologs, \%h_tmp;
    $nb ++;
}




#
#  load the exons
#

$ta->loadFile("DATA/" . $a_ref_ref->[2]);
my $a_ref = $ta->getArray();

my %CHR            = ();
my %EXONS          = ();
my %BOUNDARIES     = ();
my %BOUNDARIES_CDS = ();
my %STRAND         = ();

foreach my $r (@$a_ref) {
    my @a_exon = ($r->[2], $r->[3]);
    push @{ $EXONS{ $r->[0] } }, \@a_exon;
    $CHR{ $r->[0] }             = $r->[1]; 
    $BOUNDARIES{ $r->[0] }->[0] = (defined($BOUNDARIES{ $r->[0] }->[0])?Sets::min($BOUNDARIES{ $r->[0] }->[0], $r->[2]):$r->[2]);
    $BOUNDARIES{ $r->[0] }->[1] = (defined($BOUNDARIES{ $r->[0] }->[1])?Sets::max($BOUNDARIES{ $r->[0] }->[1], $r->[3]):$r->[3]);
    $STRAND{ $r->[0] }          = $r->[4];

    # no need to do that if this exon is not part of the CDS
    # next if (($r->[5] eq "") && ($r->[6] eq ""));
    
    if (!defined($r->[5])) {
      $r->[5] = $r->[2];
      $r->[6] = $r->[3];
    }

    
    $BOUNDARIES_CDS{ $r->[0] }->[0] = (defined($BOUNDARIES_CDS{ $r->[0] }->[0])?Sets::min($BOUNDARIES_CDS{ $r->[0] }->[0], $r->[5]):$r->[5]);
    $BOUNDARIES_CDS{ $r->[0] }->[1] = (defined($BOUNDARIES_CDS{ $r->[0] }->[1])?Sets::max($BOUNDARIES_CDS{ $r->[0] }->[1], $r->[6]):$r->[6]);
}



# 
#  get the 3'UTR
#
my $start = undef;
my $end   = undef;
my $seq   = undef;


my $utr   = $ARGV[1];



if ($STRAND{ $ARGV[0] } == 1) {

  if ($utr eq "D") {
    $start = $BOUNDARIES{ $ARGV[0] }->[1];
    $end   = $start + $len;
  } else {
    $end   = $BOUNDARIES    { $ARGV[0] }->[0]; 
    $start = $end   - $len;
  } 

  $se->setBlastDB("DATA/" . $a_ref_ref->[1]);
  $seq = $se->getSequenceFromBlastDB($CHR{ $ARGV[0] }, $start, $end);

} else {

  if ($utr eq "D") {
    $end   = $BOUNDARIES    { $ARGV[0] }->[0]; 
    $start = $end   - $len;
  } else {
    $start = $BOUNDARIES    { $ARGV[0] }->[1];
    $end   = $start + $len;
  }

  $se->setBlastDB("DATA/" . $a_ref_ref->[1]);
  $seq = $se->getSequenceFromBlastDB($CHR{ $ARGV[0] }, $start, $end);
  $seq = Sets::getComplement($seq);

}

die "pb, seq is too short\n" if (length($seq) < 10);

open OUT, ">$ARGV[0].seq";

print OUT ">$a_ref_ref->[0]\n$seq\n";

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

      if ($utr eq "D") {
	$st = $o->{ORTHOLOGS}->{ $ARGV[0] }->[ 3 ] - $offset; 
	$en = $o->{ORTHOLOGS}->{ $ARGV[0] }->[ 3 ] + $offset + $len;
      } else {
	$st = $o->{ORTHOLOGS}->{ $ARGV[0] }->[ 2 ] - $offset - $len; 
	$en = $o->{ORTHOLOGS}->{ $ARGV[0] }->[ 2 ] + $offset;
      }
      
      $st = 0 if ($st < 0);
      $seq_orth = $se->getSequenceFromBlastDB($ch, $st, $en); 

    } else {

      if ($utr eq "D") {
	$st = $o->{ORTHOLOGS}->{ $ARGV[0] }->[ 2 ] - $offset - $len;
	$en = $o->{ORTHOLOGS}->{ $ARGV[0] }->[ 2 ] + $offset;
      } else {
	$st = $o->{ORTHOLOGS}->{ $ARGV[0] }->[ 3 ] - $offset;
	$en = $o->{ORTHOLOGS}->{ $ARGV[0] }->[ 3 ] + $offset + $len;
      }
      
      $st = 0 if ($st < 0);
      $seq_orth = $se->getSequenceFromBlastDB($ch, $st, $en); 
      $seq_orth = Sets::getComplement($seq_orth);
      
    }
    if (defined($seq_orth)) {
	print OUT ">$o->{NAME}\n$seq_orth\n";
    }
    
} 
   

close OUT;


my $todo = $df->get("DIALIGN") . " -n $ARGV[0].seq";
print "$todo\n";
system($todo);


