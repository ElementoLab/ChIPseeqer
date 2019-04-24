BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my $nc    = $ta->getNbColumns();

my @COLS = ();

push @COLS, $ta->getColumn(0);



for (my $i=1; $i<$nc; $i++) {
  
  my $a_ref_col = $ta->getColumn($i);
  
  my $m  = shift @$a_ref_col;
  
  my @a_col_sorted = sort { $a <=> $b } @$a_ref_col; 
  
  my $n  = scalar(@a_col_sorted);

  my $s1 = undef;
  my $s2 = undef;
  my $std = Sets::stddev ($a_ref_col);
  my $avg = Sets::average($a_ref_col);
  if ($ARGV[1] eq "2std") {
    $s2 = $avg + 2 * $std;
  } elsif ($ARGV[1] eq "3std") {
    $s2 = $avg + 3 * $std;
  } else {
    $t = $ARGV[1];
    $s1 = $a_col_sorted[ $t      ];
    $s2 = $a_col_sorted[ $n - $t ];
  }
  
 

  my @a_col_thre   = ();
  push @a_col_thre, $m;
  foreach my $r (@$a_ref_col) {
    
    if (!defined($ARGV[2])) {

      if ($r >= $s2) {
	push @a_col_thre, 2;
      } elsif ($r <= $s1) {
	push @a_col_thre, 0;
      } else {
	push @a_col_thre, 1;
      }

    } else {
      
      if ($r >= $s2) {
	push @a_col_thre, 1;
      } else {
	push @a_col_thre, 0;
      }

    }

  }
  
  push @COLS, \@a_col_thre;
}

my $m = scalar(@COLS);
my $n = scalar(@{ $COLS[0] });
for (my $i=0; $i<$n; $i++) {
  
  print $COLS[0][$i];
  for (my $j=1; $j<$m; $j++) {
    print "\t$COLS[$j][$i]";
  
  }
  print "\n";
}
