BEGIN{ $home = `echo \$HOME`; chomp $home}

use lib "$home/PERL_MODULES";

use Sequence;

my $se = Sequence->new;
#$se->setVerbose(1);

open IN, $ARGV[0];
while (my $l = <IN>) {
    chomp $l;

    my @a = split /\t/, $l;

    my $n = $a[4];
    my $d = $a[2];

    $se->setBlastDB("$home/PROGRAMS/PHENOTYPES/GENOMES/$d/genome.faa");

    my $seq = $se->getSequenceFromBlastDB($n, 0, 0);

    if ($seq) {
        print ">$n" . " $d\n$seq\n";

    }

}
