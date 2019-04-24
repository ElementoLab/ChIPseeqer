#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Fasta;
use Getopt::Long;

my $tssfile    = undef;
my $rdm        = 'rdm';
my $genome     = undef;
my $annotation = undef;


if (@ARGV == 0) {
  die "usage: perl gene.. --tssfile --annotation --genome --rdm\n";
}

GetOptions ('tssfile=s'      => \$tssfile,
	    'annotation=s'   => \$annotation,
	    'rdm=s'          => \$rdm,
	    'genome=s'       => \$genome);


my $ta = Table->new;

#
# load real TSS info
#
$ta->loadFile($tssfile);
my $a_ref = $ta->getArray();

my %CHR_GENES = ();
foreach my $r (@$a_ref) {
  push @{ $CHR_GENES{$r->[1]} }, $r;
}

#
# load gene coord info
#
$ta->loadFile($annotation);
my $a_ref_genes = $ta->getArray();
my %CHR_GENES = ();
foreach my $r (@$a_ref_genes) {
  push @{ $CHR_GENES{$r->[1]} }, $r;
}

#
# what we need to pre-calculate
#
my %CHR_LENS   = ();

#
# load chromosome info
#
my $fa = Fasta->new;
$fa->setFile($genome);

while  (my $a_ref = $fa->nextSeq()) {
  my ($n, $s)     = @$a_ref;
  $CHR_LENS{ $n } = length( $s );
}

foreach my $r (@$a_ref) {

  $r->[0] = "$r->[0]-$rdm";

  while (1) {

    my $pos_tss = int(rand( $CHR_LENS{$r->[1]} ));
    my $ingene = 0;
    foreach my $g (@{ $CHR_GENES{ $r->[1] } }) {
      if (($pos_tss >= $g->[2]) && ($pos_tss <= $g->[3])) {
	$ingene = 1;
	last;
      }
    }

    if ($ingene == 1) {
      $r->[2] = $pos_tss;
      $r->[3] = $pos_tss;
      $r->[5] = $pos_tss;
      $r->[6] = $pos_tss;
      print join("\t", @$r) . "\n";
      last;
    }
  }

}

