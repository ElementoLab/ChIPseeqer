#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use Fasta;
use Sets;
use Blat;
use strict;
use Sequence;
use Getopt::Long;

my $oligos  = undef;
my $flank   = undef;
my $chrdir  = undef;
my $match   = undef;
my $onlytop = 0;
my $chrfile = undef;
my $outfile = undef;

if (@ARGV == 0) {
  die "Usage: perl map_oligos_to_chromosomes_using_blat.pl --oligos=FILE --chrdir=DIR --flank=1000\n";
}

GetOptions ('oligos=s'    => \$oligos,
            'chrfile=s'   => \$chrfile,
	    'outfile=s'   => \$outfile,
            'flank=s'     => \$flank,
            'onlytop=s'   => \$onlytop,
	    "match=s"     => \$match);



my $mb = Blat->new;

$mb->setVerbose(0);

$mb->setQuery($oligos);
$mb->setDatabase($chrfile);

#$mb->setOutfile();

$mb->blat($outfile);

my $h_ref = $mb->getExactMatches($outfile);

foreach my $id (%$h_ref) {
  foreach my $r (@{ $h_ref->{$id} }) {
    print "$id\t$r->{HIT_NAME}\t$r->{DFROM}\t$r->{DTO}\t$r->{DSTRAND}\n";
  }
}

#$mb->cleanup();
