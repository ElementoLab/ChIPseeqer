BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";



use Genbank;
use Data::Dumper;
use strict;

my $gb = Genbank->new;

$gb->read($ARGV[0]);

my $a_ref_feat = $gb->getFeatures();

my $ta = $gb->getTaxonomy();
my $sp = $gb->getSpecies();

$sp =~ s/\ /\_/g;

my @a_ta = split /\; /, $ta; 

my %H      = ();
my %COUNTS = ();
foreach my $f (@$a_ref_feat) {
  

  if ($f->{NAME} eq "gene") { 
    
    my $a_ref_p = $gb->getPosFromPOS( $f->{POS } );
    
    my $acc     = $gb->getAccession();

    my $st = ($a_ref_p->[2]>0?1:-1);

    print "$f->{gene}\t$acc\t$a_ref_p->[0]\t$a_ref_p->[1]\t$st\n";
    
  } 
}



