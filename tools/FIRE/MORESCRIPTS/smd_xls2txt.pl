use strict;


open IN, $ARGV[0];

my $f = $ARGV[0] . ".txt";



open OUT, ">$f" or die "cannot open $f\n";

my $cnt = 0;

my $LL_index = undef;
my $LR_index = undef;

my $name = (defined($ARGV[1])?$ARGV[1]:"Locuslink ID");


while (my $l = <IN>) {
  chomp $l;
   
  next if ($l =~ /^\!/);

  my @a = split /\t/, $l, -1;

  if ($cnt == 0) {
    my $n = scalar( @a );

    for (my $i=0; $i<$n; $i++) {
      
      if (($a[$i] eq $name) || ($a[$i] eq "Description")) {
	$LL_index = $i;
      }

      if (($a[$i] eq "Log(base2) of R/G Normalized Ratio (Mean)") || ($a[$i] eq "LOG_RAT2N_MEAN")) {
	$LR_index = $i;
      }

    }
  }

  die "$ARGV[0]: LL or LR index not found\n" if (!defined($LL_index) || !defined($LR_index));

  print OUT "$a[$LL_index]\t$a[$LR_index]\n" if (($a[$LR_index] ne "") && ($a[$LL_index] ne "EMPTY") && ($a[$LL_index] ne ""));  
  $cnt ++;
}


close IN;
close OUT;
