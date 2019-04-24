use lib qw(/home/elemento/PERL_MODULES);
use Table;

my $ta = Table->new;
$ta->loadFile("/home/elemento/DATA/DROSOPHILA/INTERACTIONS/gene_names_index.txt");
my $a_ref_dic = $ta->getArray;

my %H = ();
foreach my $r (@$a_ref_dic) {
  my $n = shift @$r;
  foreach my $s (@$r) {
    $H{ $s } = $n;
  }
  
}



#
# read in interactions
#

my $c = shift @ARGV;

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

  next if ($a[0] ne $c);

  #print "CLUSTER $a[0]\n";
  
  my @b = split /, /, $a[4];
  my $n = scalar(@b);

  my %H = ();
  for (my $i=0; $i<$n-1; $i++) {
    for (my $j=$i+1; $j<$n; $j++) {
      
      if ( defined($EDGES{ $b[$i] }{ $b[$j] }) ) {
	print $H{$b[$i]} . "\t" . $H{$b[$j]} . "\n";
	#foreach my $k (@{ $EDGES{$b[$i]}{$b[$j]} }) {
        # $H{ $k } ++;
       #}
      }

    }
  }

  foreach my $k (keys(%H)) {
    print "$H{$k}\t$k\n";
  }
}

close IN;
