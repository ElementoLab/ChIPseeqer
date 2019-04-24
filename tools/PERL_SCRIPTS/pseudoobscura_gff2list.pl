BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

#use lib qw(/home/elemento/PERL_MODULES);

use Table;
use strict;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();



my %IND = ();

foreach my $r (@$a_ref) {

  next if ($r->[1] ne ".");
  
  if ($r->[2] eq "mRNA") {
    my ($ci) = $r->[8] =~ /ID\=(.+?)\-/;
    $IND { $ci } -> [0] = $ci;
    $IND { $ci } -> [1] = $r->[0];
    $IND { $ci } -> [5] = (defined($IND { $ci } -> [5])?Sets::min($IND { $ci } -> [5], $r->[3]):$r->[3]);
    $IND { $ci } -> [6] = Sets::max($IND { $ci } -> [6], $r->[4]);
    $IND { $ci } -> [4] = ($r->[6]eq'+'?1:-1);
    
  } elsif ($r->[2] eq "CDS")  {
    my ($ci) = $r->[8] =~ /ID\=(.+?)\-/;
    $IND { $ci } -> [2] = (defined($IND { $ci } -> [2])?Sets::min($IND { $ci } -> [2], $r->[3]):$r->[3]);
    $IND { $ci } -> [3] = Sets::max($IND { $ci } -> [3], $r->[4]);
  
  }
}


#perl /home/elemento/PERL_MODULES/SCRIPTS/pseudoobscura_gff2list.pl dpse-2-r2.0.gff  | wc -l

foreach my $r (keys(%IND)) {
  print join("\t", @{ $IND{$r} }); print "\n";
}
