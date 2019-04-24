BEGIN{ $home = `echo \$HOME`; chomp $home}

use lib "$home/PERL_MODULES";

use Sets;
use Table;



my $ta = Table->new;
$ta->loadFile($ARGV[1]);
my $h_ref = $ta->getIndex(0);

open IN, $ARGV[0];
while (my $l = <IN>) {
  chomp $l;  
  my @a = split /\t/, $l, -1;
  
  if (defined($h_ref->{ $a[ $ARGV[2] ] })) {	
    print join("\t", @a);
    print "\t";
    print join("\t", @{ $h_ref->{ $a[ $ARGV[2] ] } });
    print "\n";
  } else {
    #if (defined($ARGV[3])) {
      print join("\t", @a); print "\n"; #iprint ("\tNA" x $ARGV[3]); print "\n";
    #}
  }  
}
close IN;


