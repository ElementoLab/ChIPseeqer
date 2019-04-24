use lib qw(/home/elemento/PERL_MODULES);
use strict;


my %ROWS = ();

open IN, $ARGV[0];
my $l = <IN>;
print $l;
my $cnt = 0;
while (my $l = <IN>) {
    chomp $l;
    my @a = split /\t/, $l, -1;
    my $n = shift @a;
    if (!defined($ROWS{ $n })) {
      print "$l\n";
    }
    $ROWS{ $n } = 1;
}

close IN;

