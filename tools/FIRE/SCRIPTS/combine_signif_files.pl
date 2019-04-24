use lib "$ENV{FIREDIR}/SCRIPTS";

use strict;
use Sets;
use Getopt::Long;

my $expfile=undef;
my $gap1=undef;
my $gap2=undef;
my $motiffile_dna=undef;
my $motiffile_rna=undef;

GetOptions('expfile=s'       => \$expfile,
	   'gap1=s'          => \$gap1,
	   'gap2=s'          => \$gap2,
	   'motiffile_dna=s' => \$motiffile_dna,
	   'motiffile_rna=s' => \$motiffile_rna);

my $metadir = $expfile."\_META" ;

my @motifs ;
open I, $motiffile_dna;
while(<I>)
{
    chomp;
    push(@motifs, $_);
}
close I;

open I, $motiffile_rna;
while(<I>)
{
    chomp;
    push(@motifs, $_);
}
close I;
my $fn = Sets::filename($expfile);
mkdir "$metadir/$fn\_FIRE" if (! -d "$metadir/$fn\_FIRE");
system("cp $metadir/$fn\_FIRE/DNA/$fn.signif.motifs.rep $metadir/$fn\_FIRE/DNA_RNA/$fn.signif.dna.motifs.rep");
system("cp $metadir/$fn\_FIRE/RNA/$fn.signif.motifs.rep $metadir/$fn\_FIRE/DNA_RNA/$fn.signif.rna.motifs.rep");

open OUT_dna, "> $metadir/$fn\_FIRE/DNA/$fn.signif" or die "Couldn't: $metadir/$fn\_FIRE/DNA/$fn.signif";
open OUT_dnarna_dna, "> $metadir/$fn\_FIRE/DNA_RNA/$fn.signif.dna" or die "Couldn't: $metadir/$fn\_FIRE/DNA_RNA/$fn.signif.dna";
open OUT_rna, "> $metadir/$fn\_FIRE/RNA/$fn.signif" or die "Couldn't: $metadir/$fn\_FIRE/RNA/$fn.signif";
open OUT_dnarna_rna, "> $metadir/$fn\_FIRE/DNA_RNA/$fn.signif.rna" or die "Couldn't: $metadir/$fn\_FIRE/DNA_RNA/$fn.signif.rna";

open NAME_dna, "> $metadir/$fn\_FIRE/DNA/$fn.motifnames" or die "Couldn't: $metadir/$fn\_FIRE/DNA/$fn.motifnames";
open NAME_rna, "> $metadir/$fn\_FIRE/RNA/$fn.motifnames" or die "Couldn't: $metadir/$fn\_FIRE/RNA/$fn.motifnames";

foreach my $i($gap1..$gap2)
{
    my $summaryfile = "$expfile\_g$i\_FIRE/DNA/$fn\_g$i.signif";
    open I, "< $summaryfile" or print "Couldn't: $summaryfile\n" ;
    while(<I>)
    {
	chomp ;
	my @a = split(/\t/, $_);
	next if (! (grep {$_ eq $a[0]} (@motifs)));
	print OUT_dna join("\t", @a), "\n";
	print OUT_dnarna_dna join("\t", @a), "\n";
    }
    close I;
    $summaryfile = "$expfile\_g$i\_FIRE/RNA/$fn\_g$i.signif";
    open I, "< $summaryfile" or print "Couldn't: $summaryfile\n" ;
    while(<I>)
    {
        chomp ;
        my @a = split(/\t/, $_);
        next if (! (grep {$_ eq $a[0]} (@motifs)));
        print OUT_rna join("\t", @a), "\n";
	print OUT_dnarna_rna join("\t", @a), "\n";
    }
    close I;

    my $namesfile = "$expfile\_g$i\_FIRE/DNA/$fn\_g$i.motifnames";
    open I, "< $namesfile" or print "$namesfile\n" ;
    while(<I>)
    {
        chomp ;
        my @a = split(/\t/, $_);
        next if (! (grep {$_ eq $a[0]} (@motifs)));
        print NAME_dna join("\t", @a), "\n";
    }
    close I;
    $namesfile = "$expfile\_g$i\_FIRE/RNA/$fn\_g$i.motifnames";
    open I, "< $namesfile" or print "$namesfile\n" ;
    while(<I>)
    {
        chomp ;
        my @a = split(/\t/, $_);
        next if (! (grep {$_ eq $a[0]} (@motifs)));
        print NAME_rna join("\t", @a), "\n";
    }
    close I;
}

