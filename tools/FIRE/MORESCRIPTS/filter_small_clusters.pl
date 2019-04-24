BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";
use strict;

use Table;
use Sets;
use Fasta;

my $fastafile       = Sets::get_parameter(\@ARGV, "-fastafile");
my $expfile         = Sets::get_parameter(\@ARGV, "-expfile");
my $nmin            = Sets::get_parameter(\@ARGV, "-nmin");

my $fa = Fasta->new;
$fa->setFile($fastafile);

my %H = ();
while (my $a_ref = $fa->nextSeq()) {
  my ($n, $s) = @$a_ref;
  $H{ $n } = 1;
}

my $ta = Table->new;
$ta->loadFile($expfile);
my $a_ref = $ta->getArray();

my $cnt = 0 ;
my %C = ();
foreach my $r (@$a_ref) {
  if ($cnt == 0) {
    print join("\t", @$r); print "\n";
  } else {
    
    push @{ $C{ $r->[1] } }, $r->[0] if (defined($H{ $r->[0] }));
    

    
  }

  
  $cnt ++;
}

foreach my $c (keys(%C)) {
  my $n = scalar( @{$C{$c}} );
  
  if ($n >= $nmin) {
    
    foreach my $g (@{$C{$c}}) {
      print "$g\t$c\n";
    }
    
  }
}
