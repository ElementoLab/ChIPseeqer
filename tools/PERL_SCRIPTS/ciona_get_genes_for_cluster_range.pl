use lib qw(/home/elemento/PERL_MODULES);
use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $h_ref = $ta->getIndexKV(0,1);

open IN, $ARGV[1];
<IN>; <IN>; <IN>;
while (my $l=<IN>) {
  
  my @a = split /\t/, $l, -1;

  if ($a[1] eq $ARGV[2]) {
    $start = 1;
  }
  
  if ($a[1] eq $ARGV[3]) {
    $end = 1;
  }

  if ($start == 1) {
    print "$h_ref->{$a[1]}\n" if ($h_ref->{$a[1]} !~ /NOTHING/);
  }

  if ($end   == 1) {
    exit;
  }

}
