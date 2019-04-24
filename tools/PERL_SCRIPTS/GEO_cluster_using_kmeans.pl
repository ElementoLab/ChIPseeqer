#!/usr/bin/perl

print "\nDetermine # clusters ... \n";
$todo = "getnbclusters $ARGV[0]";
my $nbc = `$todo`; chomp $nbc;

print "$nbc.\n";


print "Run cluster 3.0 ... ";

$todo = "/usr/local/bin/cluster -cg a -ng -g 2 -k $nbc -r 10 -f $ARGV[0]";
system($todo) == 0 or die "Could not cluster ... \n";

print "Done (clustered into $nbc clustered)\n";
