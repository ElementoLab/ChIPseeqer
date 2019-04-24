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
  foreach my $fh (@FH) {
    my $l = <$fh>;
    if (!defined($l)) {
      $o = 1;
      last;
    }

    chomp $l;
    my @b = split /\t/, $l, -1;
    if ($cnt == 0) {
      print "$b[0]";
    }
    
    if ($lcnt != 0) {
      print "\t$b[1]";
    } else {
      print "\t$ARGV[$cnt]";
    }
    $cnt ++;
  }
  $lcnt ++;
  if ($o == 1) {
    last;
  }
  print "\n";
}
  



