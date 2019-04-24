BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;

my $a_ref_cols = Sets::readSet($ARGV[0]);



open IN, $ARGV[1] or die "Cannot open $ARGV[1]\n";

my $l = <IN>;
chomp $l;

my @order = ();
my @h = split /\t/, $l, -1;
my @hc = @h;
shift(@h);

# build an order in which to visit the columns
push @order, 0;


my $i = 1;
my %H = ();
foreach my $r (@h) {
  $H{ $r } = $i;
  $i ++;
}

foreach my $r (@$a_ref_cols) {
  push @order, $H{$r};
}


my @newa = ();
for (my $i=0; $i<@order; $i++) {
  push @newa, $hc[$order[$i]];
}

print join("\t", @newa) . "\n";

while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;

  my @newa = ();
  for (my $i=0; $i<@order; $i++) {
    push @newa, $a[$order[$i]];
  }

  print join("\t", @newa) . "\n";
}

close IN;
