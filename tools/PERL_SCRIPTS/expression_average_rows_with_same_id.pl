#!/usr/bin/perl



use strict;


my %ROWS = ();

open IN, $ARGV[0];
my $l = <IN>;
print $l;
my $cnt = 0;
while (my $l = <IN>) {
    chomp $l;
    my @a = split /\t/, $l, -1;
    
    my $n = shift @a;
    
    push @{ $ROWS{ $n } }, \@a;

    $cnt ++;

    #last if ($cnt == 100);
}

close IN;


my @LL = sort(keys(%ROWS));

foreach my $k (@LL) {
    
    my $n = scalar(@{ $ROWS{ $k } });
    my $l = scalar(@{ $ROWS{ $k }->[0] });
    
    #print "n=$n\n";

    my @exp = ();
    for (my $i=0; $i<$l; $i++) {
      my $cnt = 0;
      my $sum = undef;
      for (my $j=0; $j<$n; $j++) {
	if (($ROWS{ $k }->[$j]->[$i] ne "") && ($ROWS{ $k }->[$j]->[$i] ne 'NULL')) {
	  $sum += $ROWS{ $k }->[$j]->[$i];
	  $cnt += 1;
	}
      }
      if ($cnt > 0) {
	my $prec = 1;
	if ($ENV{PREC} ne "") {
	  $prec = $ENV{PREC};
	}
	$exp[$i] = sprintf("%3.$prec" . "f", $sum / $cnt);
      } else {
	$exp[$i] = "";
      }
      
      
    }
    
    print "$k\t"; print join("\t", @exp); print "\n";

}
