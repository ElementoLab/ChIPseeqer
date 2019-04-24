BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;
use Table;
use strict;



open IN, $ARGV[0];
my $kmer = $ARGV[1];
srand;

if ($ARGV[2] == 1) {
  my $a = Sets::shuffle_re($kmer);
  $kmer = join("", @$a);

  print "$kmer\n";
}

my $h_ref = undef;
if (defined($ARGV[3])) {
  my $ta = Table->new;
  $ta->loadFile($ARGV[3]);
  $h_ref = $ta->getIndexKV(0,1);
}

my $orf  = undef;

my $cnt_ok = 0;
my $cnt_th = 0;
while (my $l = <IN>) {
  chomp $l;
  
  next if ($l eq "");

  #print "$l\n";

  if ($l =~ /^\>(.+)$/) {

    $orf = $1;
    
    my $lh = <IN>;
 
    my @o  = ();
    my $lm = <IN>; chomp $lm; push @o, $lm;
    my $lr = <IN>; chomp $lr; push @o, $lr;
    my $lc = <IN>; chomp $lc; push @o, $lc;
    
    if (defined($h_ref) && ($h_ref->{$orf} == 0)) {
      next;
    }

    my $kmer_c = Sets::getComplement($kmer);

    while ($lh =~ /($kmer|$kmer_c)/ig) {
      
      my $r = $1;
      my $p = pos($lh) - length($&);
      $cnt_th ++ ;

      #print "$p\n";

      my $l = length($r);

      my $cnt = 0;
      foreach my $s (@o) {
	my $ss =  substr($s, $p, $l); 
	if ($ss =~ /($kmer|$kmer_c)/i) {
	  $cnt ++;
	}
      }
      if ($cnt == 3) {
	$cnt_ok ++ ;
      }

    }
  }

}
close IN;

my $r =  sprintf("%5.4f", $cnt_ok/$cnt_th);
print "$cnt_ok/$cnt_th\t$r\n";
