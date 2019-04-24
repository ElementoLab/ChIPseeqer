use Getopt::Long;

if (@ARGV == 0) {
  die "Args --matches=FILE --pattern=STR --col=INT\n";
}

my $matches = undef;
my $pattern = undef;
my $col     = 0;

GetOptions("matches=s" => \$matches,
	   "col=s"     => \$col,
           "pattern=s" => \$pattern);

open IN, $matches or die "Cannot open $matches\n";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  if ($a[$col] =~ /$pattern/) {
    print "$l\t1\n";
  } else {
    print "$l\t0\n";
  }
}
close IN;

