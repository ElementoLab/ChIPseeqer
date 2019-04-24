BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Genbank;
use Data::Dumper;
use strict;



my %CODON_TABLE = (
		   TCA => 'S',TCG => 'S',TCC => 'S',TCT => 'S',
		   TTT => 'F',TTC => 'F',TTA => 'L',TTG => 'L',
		   TAT => 'Y',TAC => 'Y',TAA => '*',TAG => '*',
		   TGT => 'C',TGC => 'C',TGA => '*',TGG => 'W',
		   CTA => 'L',CTG => 'L',CTC => 'L',CTT => 'L',
		   CCA => 'P',CCG => 'P',CCC => 'P',CCT => 'P',
		   CAT => 'H',CAC => 'H',CAA => 'Q',CAG => 'Q',
		   CGA => 'R',CGG => 'R',CGC => 'R',CGT => 'R',
		   ATT => 'I',ATC => 'I',ATA => 'I',ATG => 'M',
		   ACA => 'T',ACG => 'T',ACC => 'T',ACT => 'T',
		   AAT => 'N',AAC => 'N',AAA => 'K',AAG => 'K',
		   AGT => 'S',AGC => 'S',AGA => 'R',AGG => 'R',
		   GTA => 'V',GTG => 'V',GTC => 'V',GTT => 'V',
		   GCA => 'A',GCG => 'A',GCC => 'A',GCT => 'A',
		   GAT => 'D',GAC => 'D',GAA => 'E',GAG => 'E',
		   GGA => 'G',GGG => 'G',GGC => 'G',GGT => 'G');


my $gb = Genbank->new;

$gb->read($ARGV[0]);

my $a_ref_feat = $gb->getFeatures();

my $ta = $gb->getTaxonomy();
my $sp = $gb->getSpecies();

print "$sp\t/\t$ta\n";

$sp =~ s/\ /\_/g;

my @a_ta = split /\; /, $ta; 

my $file = $ARGV[0]; $file =~ s/\.txt//;
$file = $file . "_$a_ta[1].fa"; 
die "NO\n" if ($file eq $ARGV[0]);

#open OUT, ">$file";  

#print "$sp\t$ta\n";

my %H      = ();
my %COUNTS = ();


foreach my $f (@$a_ref_feat) {
  
  if ($f->{NAME} eq "CDS") { 
    
    #print "$f->{POS}\n";
    my $seqnt = $gb->getSubseqFromPOS( $f->{POS } );
    my $seqaa = $f->{translation}; 
    my $l     = length($seqnt);
    
    my @a = split //, $seqnt;
    my @b = split //, $seqaa;

    for (my $i=0, my $j=0; $i<$l; $i+=3, $j++) {
      my $codon = $a[$i] . $a[$i+1] . $a[$i+2];
      my $aa    = $b[$j];
      if (!$aa) {
	$aa = '?';
      }

      $COUNTS{ $aa }{ $codon } ++;
      
    }
  }
  
}

#foreach my $aa (keys(%COUNTS)) {
my $aa = 'W';
foreach my $co (keys(%{ $COUNTS{$aa} })) {
  my $t = Sets::translate($co);
  print "$aa\t$co\t$COUNTS{$aa}{$co}\t$t\n";
}
#}
