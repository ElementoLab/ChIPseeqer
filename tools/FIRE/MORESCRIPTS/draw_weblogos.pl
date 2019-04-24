BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Sets;

my $f = $ARGV[0];
my $d = $ARGV[1]; die "Please specify d/\n" if (!defined($d)); 
mkdir $d if (! -e $d);
my $rootdir = ".";

my $cnt = 0;

open IN, $f;
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l;
  
  my $mo = Sets::myre2wm($a[0]);
  open OUT, ">$d/$cnt.txt" or die "cannot open $d/$cnt.txt\n";
  print OUT $mo;
  close OUT;
  
  system("$rootdir/weblogo/seqlogo -f $d/$cnt.txt -F EPS  -a -c -M -n -Y -w 5 -h 3 > $d/$cnt.eps");
  $cnt ++;
}
close IN;
