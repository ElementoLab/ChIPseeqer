BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Fasta;
use strict;

my $ta = Table->new;


# load exp file
my $expfile  = $ARGV[0];
$ta->loadFile($expfile);
my $a_ref = $ta->getArray();

my $fastafile = $ARGV[1];


my $fa = Fasta->new;
$fa->setFile($fastafile);

my %SEQ = ();
while (my $a_ref = $fa->nextSeq()) {
    my ($n, $s) = @$a_ref;
    $SEQ{$n} = $s;
}

my $d = "$expfile\_LEAVEONEOUT";
mkdir $d if (! -e $d);

my $cnt = 0;
shift @$a_ref;
foreach my $r (@$a_ref) {
  
  open OUT, ">$d/$cnt.txt";
  foreach my $s (@$a_ref) {
    if ($r->[0] ne $s->[0]) {
      print OUT join("\t", @$s) . "\n";
    } else {
      print "$r->[0]\t";
    }
  }
  close OUT;
  
  my $todo = "perl fire.pl --maxselected=1 --shuffle_mifind=100 --expfiles=$d/$cnt.txt --exptype=discrete --fastafile_dna=BICOID/124_Dmel_Enc.fa  --maxfreq=0.5 --gap=2 --k=8 --optimslow=1 --dodef=0 --doremovedups=1 --domifind=1 --domioptimize=1 --dorna=0 --dodnarna=0 --submit=1 > /dev/null";
 
  #my $todo = "perl fire.pl --expfiles=$d/$cnt.txt --exptype=discrete --fastafile_dna=BICOID/124_Dmel_Enc.fa  --maxfreq=0.5 --dorna=0 --dodnarna=0 --dodef=0 --doremovedups=1 --domifind=1 --domioptimize=1 > /dev/null";
 
 
  system($todo);


  my $in = "$d/$cnt.txt\_FIRE/DNA/$cnt.txt.optim";
  open IN, $in;
  my $l = <IN>;
  my @a = split /\t/, $l;
  
  my $re = $a[0];
  print "$re\t";

  my $answer = $r->[1];

  my $sc = Sets::getComplement($SEQ{$r->[0]});
  my $guess = 0;
  if (($SEQ{$r->[0]} =~ /$re/) || ($sc =~ /$re/)) {
    $guess = 1;
  }

  print "$guess\t$answer\t";


  if ($guess == $answer) {
    print "WIN\t";
    if ($guess == 0) {
      print "TN\n";
    } else {
      print "TP\n";
    }	
    
  } else {
    print "LOSE\t";
    if (($guess == 1) && ($answer == 0)) {
      print "FP\n";
    } else {
      print "FN\n";
    }
  }

  $cnt ++;
}
  
