use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my $r = shift @$a_ref;
print join("\t", @$r) . "\n";
my %H = ();
foreach my $r (@$a_ref) {
  $H{$r->[0]} = $r->[1];
}

foreach my $k (keys(%H)) {
  print "$k\t$H{$k}\n";
}

