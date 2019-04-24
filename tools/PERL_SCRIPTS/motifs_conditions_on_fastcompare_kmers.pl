my $home;
BEGIN{ $home = `echo \$HOME`; chomp $home}

use lib "$home/PERL_MODULES";
use Table;
use Sets;
use Getopt::Long;
use strict;

my $kmerfile = undef;
my $expfile  = undef;
my $seqfile  = undef;
my $bonfplus = undef;

GetOptions ('kmerfile=s'   => \$kmerfile,
	    'expfile=s'    => \$expfile,
	    'seqfile=s'    => \$seqfile,
	    'bonfplus=s'   => \$bonfplus);


die "$expfile does not exist\n" if (! -e $expfile);

my $ta = Table->new;
$ta->loadFile($kmerfile);
my $h_ref_kmers = $ta->getIndex(0);
my $a_ref_kmers = $ta->getColumn(0);
my $n           = scalar(@$a_ref_kmers);
my $b           = $n * $bonfplus;
my %H           = ();

foreach my $kmer (@$a_ref_kmers) {
    if ($kmer =~ /^([ATCG]+)(N+)([ATCG]+)$/) {
	push @{ $H{ length($1) . "_" . length($2) . "_" . length($3) } }, "$1$3";
    } else {
	push @{ $H{ length($kmer) } }, $kmer;
    }
}

foreach my $k (sort(keys(%H))) {
    print "k=$k\n";
    my $file = "/tmp/toto";
    Sets::printSet(\@{ $H{$k} } );
    Sets::writeSet(\@{ $H{$k} }, $file);
    my $myk = undef;
    if ($k =~ /^(\d)\_(\d)\_(\d)$/) {
	$myk = $1 + $3;
    } else {
	$myk = $k;
    }

    my $todo = "$home/PROGRAMS/MOTIFS_CONDITIONS/motifs_conditions_smallset -expfile $expfile -fastafile $seqfile -kmersize $myk -regul up -singlestrand 1 -kmerfile $file -bonfplus $b";
    if ($k =~ /^(\d)\_(\d)\_(\d)$/) {
	$todo .= " -k1 $1 -gap $2 ";
    }
    $todo .= " | tee file.$k.up.txt";
    print "$todo\n";
    system("$todo") == 0 or die "cannot execute 1 \n";
    
    open IN, "file.$k.up.txt";
    while (my $l = <IN>) {
	chomp $l;
	if ($l !~ /Cond/) {
	    my @a = split /\t/, $l;
	    $h_ref_kmers->{ $a[0] }->[5] ++;
	}
    }
    close IN;

    my $todo = "$home/PROGRAMS/MOTIFS_CONDITIONS/motifs_conditions_smallset -expfile $expfile -fastafile $seqfile -kmersize $myk -regul down -singlestrand 1 -kmerfile $file -bonfplus $b";
    if ($k =~ /^(\d)\_(\d)\_(\d)$/) {
	$todo .= " -k1 $1 -gap $2 ";
    }
    $todo .= " | tee file.$k.do.txt";    
    system("$todo") == 0 or die "cannot execute 2 \n";
    
    open IN, "file.$k.do.txt";
    while (my $l = <IN>) {
	chomp $l;
	if ($l !~ /Cond/) {
	    my @a = split /\t/, $l;
	    $h_ref_kmers->{ $a[0] }->[6] ++;
	}
    }
    close IN;

}

my $u = 0;
my $d = 0;
foreach my $kmer (@$a_ref_kmers) {
    $u ++ if ($h_ref_kmers->{ $kmer }->[5]);
    #$h_ref_kmers->{ $kmer }->[5] = 0 if (!defined($h_ref_kmers->{ $kmer }->[5]));

    #$h_ref_kmers->{ $kmer }->[6] = 0 if (!defined($h_ref_kmers->{ $kmer }->[6]));
    $d ++ if ($h_ref_kmers->{ $kmer }->[6]);
    print sprintf("%11s\t", $kmer);
    shift @{ $h_ref_kmers->{ $kmer } };
    print join("\t", @{ $h_ref_kmers->{ $kmer } }); print "\n";

}

print "$u\t$d\n";
