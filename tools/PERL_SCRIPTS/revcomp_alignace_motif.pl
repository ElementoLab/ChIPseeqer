use lib "$ENV{HOME}/PERL_MODULES";
use Sets;

open IN, $ARGV[0];

while (my $l = <IN>) {
 chomp $l;
 if ($l !~ /[^ACGT]/) {
   $l = Sets::getComplement($l);
 }
 print "$l\n";
}

close IN;
