BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use Table;
use Fasta;
use Sets;
use strict;

# 1. load ENSG to save, if they have a LLID. 
my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my %H = ();
my $a_ref = $ta->getArray();

foreach my $r (@$a_ref) {
  if ($r->[1] ne "") {
    push @{ $H{ $r->[0] } }, $r->[1] if (!Sets::in_array($r->[1], @{ $H{ $r->[0] } }));
    #print "pushing $r->[1] on $r->[0]\n";
  }
}



my $fa = Fasta->new;
$fa->setFile($ARGV[1]);

while (my $a_ref = $fa->nextSeq()) {
    my ($n, $s) = @$a_ref;
    if (defined($H{$n})) {
      
      foreach my $v (@{ $H{ $n } }) {
	print ">$v\n$s\n\n";
      }
      
    }
}
