BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use FileHandle;

my @FH = ();
foreach my $f (@ARGV) {
  my $fh = IO::File->new($f);
  push @FH, $fh;
}


my $lcnt = 0;
while (1) {
  my @a = ();
  my $o = 0;
  my $cnt = 0;
  
  my @G = ();
  my @P = ();
  my $h = undef;
  foreach my $fh (@FH) {
    my $l = <$fh>;
    if (!defined($l)) {
      $o = 1;
      last;
    }

    chomp $l;
    my @b = split /\t/, $l, -1;
    if ($cnt == 0) {
      $h = $b[0];
    }
    
    if ($lcnt != 0) {
      push @G, $b[1];
      push @P, $b[2] if ($b[2] eq 'P');
    } else {
      print "\t$ARGV[$cnt]";
    }
    $cnt ++;
  }
  $lcnt ++;
  if ($o == 1) {
    last;
  }

  if (scalar(@P) > (scalar(@G)/2.0)) {
    print "$h\t" . join("\t", @G) . "\n";
  }

}
  



