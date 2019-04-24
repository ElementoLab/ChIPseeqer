
use lib "$ENV{FIREDIR}/SCRIPTS";
use strict;
use Sets;
use Getopt::Long;


if (@ARGV==0)
{
	print "perl combine_gapped_motifs.pl --metadir=DIR --expfile=FILE --gap=X-Y\n";
	die;
}

my $expfile=undef;
my $gap=undef;
my $metadir=undef;
GetOptions ('metadir=s'	             => \$metadir,
            'expfile=s'              => \$expfile,
	    'gap=s'                  => \$gap);

my @motifs_rna ;
my @motifs_dna ;
my $gap1 = $gap; my $gap2=$gap ;
if ($gap =~ /-/)
{
    $gap =~ /(\S+)-(\S+)/ ;
    $gap1 = $1 ;
    $gap2 = $2 ;
}
my @motifs_dna ;
my @motifs_rna ;
foreach my $g ($gap1..$gap2)
{
    my $dir = "$expfile\_g$g"."\_FIRE" ;
    my $fn = Sets::filename($expfile) ;
    $fn .= "_g$g";
    my $f = $dir."/DNA/$fn.summary" ;
    open I, "< $f" or next ;
    while (my $l = <I>)
    {
        chomp $l ;
        my ($motif, @a) = split(/\t/, $l) ;
        push(@motifs_dna, $motif) ;
    }
    close I ;
    my $f = $dir."/RNA/$fn.summary" ;
    open I, "< $f" or next ;
    while (my $l = <I>)
    {
        chomp $l ;
        my ($motif, @a) = split(/\t/, $l) ;
        push(@motifs_rna, $motif) ;
    }
    close I ;
}
open O, "> $metadir/motifs_dna_redundant" ;
print O join("\n", @motifs_dna) ;
open O, "> $metadir/motifs_rna_redundant" ;
print O join("\n", @motifs_rna) ;
close O ;
