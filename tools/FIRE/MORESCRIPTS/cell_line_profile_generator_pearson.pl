use strict ;
use Data::Dumper ;
use Sets ;

open(DATA, "< TEST/NCI-60_parsed_averaged_varnorm.txt.2") or die;
my $L = <DATA> ;
$L =~ s/\s+$// ;
my ($dummy, @cell_lines) = split(/\t/, $L) ;
my %Data ;
my %Data_h ;
my @genes ;
my $gene_count =0 ;
while(<DATA>)
{
    s/\s+$// ;
    my($gene, @data) = split(/\t/, $_) ;
    for (my $i=0 ; $i<@data ; $i++)
    {
	push(@{$Data{$cell_lines[$i]}}, $data[$i]) ;
    }
    
    push(@genes, $gene) ;
    $gene_count++ ;
}

open(CON, "< TEST/doubling_time.txt") or die;
<CON> ;
my %cond_profile ;
while(<CON>)
{
    s/\s+$// ;
    my ($cell, @a) = split(/\t/, $_) ;
    $cell =~ s/\s+$// ;
    next if (! defined (@{$Data{$cell}}) or $a[0] =~ /NA/) ;
    #if ($a[-3] eq "GCB")
    #{
#	$cond_profile{$cell} = 1 ;
#    }
#    else
#    {
#	$cond_profile{$cell} = 0 ;
#    }
    
    $cond_profile{$cell} = $a[0] ;

#    $cond_profile{$cell} = 1 if $a[15] eq "I" ;
#    $cond_profile{$cell} = 2 if $a[15] eq "II" ;
#    $cond_profile{$cell} = 3 if $a[15] eq "III" ;
#    $cond_profile{$cell} = 4 if $a[15] eq "IV" ;

    print $cell, "->", $cond_profile{$cell}, "\n" ;
}
print "press enter\n" ;
<STDIN> ;
open(OUT, "> dbl.txt") or die;
print OUT "gene\tR\n" ;
my @R ;
for (my $i ; $i<$gene_count ; $i++)
{
    my @r1 ;
    my @r2 ;
    foreach my $cell (keys %cond_profile)
    {
	push (@r1, $cond_profile{$cell}) ;
	push (@r2, $Data{$cell}->[$i]) ;
    }
    
    $R[$i] = Sets::pearson(\@r1,\@r2) ;
    print $genes[$i], "\t", $R[$i], "\n" ;
    print OUT $genes[$i], "\t", $R[$i], "\n" ;
}
