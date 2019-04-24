
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";



use Sets;
use Table;
use Fasta;
use Getopt::Long;
use strict;

my $micrornas = undef;
my $kmerfile  = undef;
my $col       = 0;
my $limit     = undef;
my $only5     = undef;

if (!$ARGV[1]) {
  die "Usage : perl match_micrornas_to_kmers.pl --micrornas=FILE --kmerfile=FILE --limit=INT --col=INT --only5=INT\n";
}



GetOptions ('micrornas=s'       => \$micrornas,
	    'kmerfile=s'        => \$kmerfile,
	    'limit=s'           => \$limit,
	    'only5=s'           => \$only5,
	    'col=s'             => \$col);


my $ta = Table->new;
if (defined($limit)) {
  $ta->setLimit($limit);
}
$ta->loadFile($kmerfile);

my $a_ref = $ta->getArray();

my $fulldisplay = 1;

my $i = 1;
foreach my $r (@$a_ref) {
    
    
  my $fa = Fasta->new;
  $fa->setFile($micrornas);

  my $match = 0;
  my $rank  = 1;
  while (my $a_seq = $fa->nextSeq()) {
    
    my ($n, $s) = @$a_seq;
    
    $s =~ s/t/u/g;
    
    my $ss = uc($s);
    $ss =~ s/U/T/g;
    
    $ss = Sets::getComplement($ss);

    my $j = 0;

    my $kmer = substr($r->[$col], 0, length($r->[$col])-$j);
    my $km   = $kmer;
	
    $km =~ s/N/\./g;


    if ($ss =~ /$km/) {
	    
      $km = Sets::getComplement($km);
      $km =~ s/N/\./g;
      
      my $lkm = lc($km);
      $lkm =~ s/t/u/g;
	    
      $s =~ s/$lkm/$km/g;

	    
      my $a_ref_pos = Sets::getREMotifPositions($km, $s);
	    
	    
      my $p = shift @$a_ref_pos;
	    
	    
      if (defined($only5) && ($p > 1)) {
	$rank++;
	next;
      }
	    
      $s =~ s/T/U/g;
      print "$i: $r->[$col]\t$rank:$n\t$s\t$p\n"; 

      $match++ if ($p > 1);

      
    }
	    
  }

  $fa->dispose();
    
  $i++;
}
