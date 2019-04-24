BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Fasta;
use Sets;
use strict;

srand(100); # 10 and 100

my $fa      = Fasta->new;
$fa->setFile($ARGV[0]);

my $a_ref          = $fa->nextSeq();
my ($n, $s)        = @$a_ref;
my ($le, $se, $re) = $s =~ /^(.)(.{6})(.+)$/;

my @ra = ($le, split(//, $re));

my $a_sh = Sets::shuffle_array(\@ra);

my $le = shift @$a_sh;

my $re = join('', @$a_sh);


print ">$n\n$le$se$re\n";
