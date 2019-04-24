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

my $file = $ARGV[0]; $file =~ s/\.txt//;
$file = $file . "_$a_ta[1].fa"; 
die "NO\n" if ($file eq $ARGV[0]);

open OUT, ">$file";  

#print "$sp\t$ta\n";

my %H      = ();
my %COUNTS = ();
foreach my $f (@$a_ref_feat) {
   
  
  
  

  if ($f->{NAME} eq "CDS") { 
    
    #print "$f->{POS}\n";

    my $seqnt =  $gb->getSubseqFromPOS( $f->{POS } );
    my $seqaa = Sets::translate($seqnt);
    $seqaa =~ s/\*$//g;
    
    


    #print "aa1=$seqaa\n";
    
    my $ff = Sets::getSequencesIdentity($seqaa, $f->{translation});
    $ff = sprintf("%4.1f", 100 * $ff);
    my $w = undef;
    $w = "***" if ($ff < 100);


    #print Dumper($f);

    my $name = undef;
    if (defined($f->{gene})) {
      $name = $f->{gene};
    } elsif (defined($f->{product})) {
      $name = $f->{product};
    }

    $name =~ s/\ /\_/g;
    
    my $tname = $name;

    if (defined($H{$name})) {
      $tname .= "_$COUNTS{$name}";     
    }

    $H     {$name} = 1;
    $COUNTS{$name} ++;

    if (defined( $f->{translation} )) {
      
      print OUT ">$tname" . "_$sp\n";
      print OUT "$f->{translation}\n";
      print OUT "\n";
    
    }
    
    #print "$ff% $w\n";

  } 
}


close OUT;
