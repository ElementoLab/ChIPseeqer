use lib qw(/home/olly/PERL_MODULES);

use Sets;

use File::Copy;
use Fasta;


open IN, $ARGV[0];
my $s = <IN>;
chomp $s;
close IN;





my @a = split /\-+/, $s;

my $cnt = 0;
foreach my $c (@a) {
    if ((length($c) >= 6) && ($c =~ /[GCgc]/)) {
	$cnt++;
    } 
} 

if ($cnt > 0) {
    print "$ARGV[0]: YES\n";
   copy($ARGV[0], "CONSERVED/$ARGV[0]");
} else {
    print "$ARGV[0] : NO\n";
} 
