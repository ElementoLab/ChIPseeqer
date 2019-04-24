use lib qw(/home/elemento/PERL_MODULES);


#
# read in interactions
#

my $s = shift @ARGV;

my %EDGES = ();

while (my $f = shift @ARGV) {
  open IN, $f;
  while (my $l = <IN>) {
    chomp $l;
    my @a = split /\t/, $l, -1;
    
    push @{ $EDGES{ $a[0] }{ $a[1] } }, $f;
    push @{ $EDGES{ $a[1] }{ $a[0] } }, $f;
    
    
  }
  close IN;
}


open IN, $s;

my $l = <IN>;
my $l = <IN>;
my $l = <IN>;
my $l = <IN>;
my $l = <IN>;
my $l = <IN>;
my $l = <IN>;
my $l = <IN>;


while (my $l = <IN>) {
    
  chomp $l;
  my @a = split /\t/, $l;
  print "CLUSTER $a[0]\n";
  
  my @b = split /, /, $a[4];
  my $n = scalar(@b);

  my %H = ();
  for (my $i=0; $i<$n-1; $i++) {
    for (my $j=$i+1; $j<$n; $j++) {
      
      if ( defined($EDGES{ $b[$i] }{ $b[$j] }) ) {
       foreach my $k (@{ $EDGES{$b[$i]}{$b[$j]} }) {
	 $H{ $k } ++;
       }
      }

    }
  }

  foreach my $k (keys(%H)) {
    print "$H{$k}\t$k\n";
  }
}

close IN;
