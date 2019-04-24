#!/usr/bin/perl

BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use Sets;
use Statistics::Test::WilcoxonRankSum;

my $ta = Table->new;

$ta->loadFile($ARGV[1]);

#
# load expression value
#
my $h_ref = $ta->getIndexKV(0,1);

#
# load 
#
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();
shift @$a_ref;

my %H = ();
foreach my $r (@$a_ref) {
  if (defined($h_ref->{ $r->[0] })) {
    push @{ $H{ $r->[1] } }, $h_ref->{ $r->[0] };
  } else {
    #print STDERR "Attention, $r->[0] in expfile.\n";
  }
}

my $wc =  Statistics::Test::WilcoxonRankSum->new( { exact_upto => 30 } );

$wc->load_data($H{"0"}, $H{"1"});

my $prob = $wc->probability();

my $pf = sprintf '%f', Sets::log10($prob); # prints 0.091022
#print "$pf\n";
my $pf = sprintf("%3e", $prob);

#print $wc->probability_status();

# prints something like:
# Probability:   0.002797, exact
# or
# Probability:   0.511020, normal approx w. mean: 104.000000, std deviation:  41.840969, z:   0.657251

#my $pstatus = $wc->probability_status();
#print "$pstatus\n";
# $pstatus is like the strings above


#$wc->probability_summary();

my $a = Sets::average($H{"0"});
my $b = Sets::average($H{"1"});

my $fe = $ARGV[1];
$fe =~ s/\.txt.*$//;
print "$fe\t";
print Sets::unTracey($ARGV[0]);
print "\t";

#if ($b > $a) {
#  print -1 * $pf;
#} else {
  print 1 * $pf;
#}

print "$a\t$b\n";
