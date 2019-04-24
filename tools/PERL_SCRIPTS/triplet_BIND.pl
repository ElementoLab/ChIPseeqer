use lib qw(/home/olly/PERL_MODULES /home/olly/PERL_MODULES/Hypergeom /home/olly/PERL_MODULES/Hypergeom/blib/lib /home/olly/PERL_MODULES/Hypergeom/blib/arch/auto/Hypergeom);
use GO_func;
use Sets;

#
#  load inter
#
open IN, $ARGV[1];
while (my $l = <IN>) {
    chomp $l;
     my @a = split /\t/, $l;
    $INTER{ $a[0] } { $a[1] } = 1;
    $INTER{ $a[1] } { $a[0] } = 1;
}
close IN;

open IN, $ARGV[0];
my $i = 0;
while (my $l = <IN>) {
    chomp $l;
    
    my @a = split /\t/, $l;

    my $r1 = $INTER{ $a[0] } { $a[1] };
    my $r2 = $INTER{ $a[1] } { $a[2] };
    my $r3 = $INTER{ $a[0] } { $a[2] };

    my $s  = $r1 + $r2 + $r3;

    if ($s == 2) {
	    print "$l\n";
    }

    $i ++;


    last if ($i == 100000);

}
close IN;
