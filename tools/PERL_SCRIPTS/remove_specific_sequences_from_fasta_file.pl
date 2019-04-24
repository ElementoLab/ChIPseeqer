use lib qw(/home/olly/PERL_MODULES);

use Fasta;
use Table;
use strict;

my $ta = Table->new;
$ta->loadFile($ARGV[1]);

my $h_ref = $ta->getIndex(0);


my $f = Fasta->new;
$f->setFile($ARGV[0]);




while ( my $a_ref = $f->nextSeq ) {
    my ($name, $seq) = @{$a_ref}; 
    my $nn = $name;
    #$nn =~ s/\.1$//;
    print ">$name\n$seq\n\n" if (!defined($h_ref->{$nn}));
}
