use lib qw(/home/elemento/PERL_MODULES);
use Table;

my $ta = Table->new;

my $h_ref = $ta->getBimensionalHash($ARGV[0]);

#
#  read CDT file now .. 
#


open IN, $ARGV[1] or die "can't open $ARGV[1] ..\n";
my $l = <IN>; print $l;
chomp $l;
my @a = split /\t/, $l, -1;
shift @a;
shift @a;
shift @a;
shift @a;

$l = <IN>; print $l; $l = <IN>; print $l;

while (my $l = <IN> ) {
  
  chomp $l;
  my @b = split /\t/, $l, -1;
    
  print shift @b; print "\t";
  print shift @b; print "\t";
  my $n = shift @b; print "$n\t";
  print shift @b;

  my $t = scalar(@b);
  for (my $i=0; $i<$t; $i++) {
    
    
    print "\t" . $h_ref->{ $n } ->{ $a[$i] }; 
    
  }
  print "\n";
}
