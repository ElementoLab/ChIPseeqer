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

    my $st = undef;

    my $st1 = undef;
    if ($f->{locus_tag} =~ /c$/) {
      $st1 = -1;
    } else {
      $st1 = 1;
    }

    my $st2 = ($a_ref_p->[2]>0?1:-1);

    my $st  = undef;
    if ($st1 == $st2) {
      $st = $st1;
    } elsif ($st1 == -1) {
      $st = -1;
    } else {
      $st = $st2;
    }

    print "$f->{locus_tag}\t$acc\t$a_ref_p->[0]\t$a_ref_p->[1]\t$st2\n";
    
    #my $seqnt = $gb->getSubseqFromPOS( $f->{POS } );
    #print "$seqnt\n";

  } 
}



