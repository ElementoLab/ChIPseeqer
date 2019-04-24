BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use Table;
use strict;

my $logfile = $ARGV[0];
my $kmer    = "TCCGTAC";
my $nr      = 1;

my $cons = $ARGV[1];


my $h_ref = undef;
if (defined($cons)) {
  my $ta = Table->new;
  $ta->loadFile($ARGV[1]);
  $h_ref = $ta->getIndexKV(0,1);
}


my $h_ref_ace = undef;
my $ace       = $ARGV[2];
if (defined($ace)) {
  my $ta = Table->new;
  $ta->loadFile($ace);
  $h_ref_ace = $ta->getIndexKV(0,1);
}




my $c = "c.txt";
open OUT, ">$c";

open IN, $ARGV[0];

my $cnt = 0;
my $mi  = undef;
my $re  = undef; 
my $oldre = undef;
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;

  

  if (($a[0] eq $kmer) && (($a[1] == $nr) || ($a[1] == -1))) {
    
    if ($a[4] ne "") {
      $mi = $a[4];
    } else {
      #
    } 
    

    if ($a[3] ne "") {
      $re = $a[3]; $oldre = $re;
      print OUT "$re\n";
      $re =~ s/\./N/g;
    } else {
      $re = "";
    }
    print "$cnt\t$mi\t$re";
    
    if (defined($cons)) {
      print "\t" . $h_ref->{$oldre}/100;
    } else {
      print "\t";
    }
    
    if (defined($ace)) {
      print "\t" . $h_ref_ace->{$oldre};
    } else {
      print "\t";
    }

    
    $cnt ++;
    print "\n";


  }
  
}

close OUT;
