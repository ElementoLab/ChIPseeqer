#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Getopt::Long;
use Sets;
use strict;

if (@ARGV == 0) {
  die "Args: --kggfile=FILE --datafile=FILE --groups=TXT\n";
}

my $kggfile  = undef;
my $datafile = undef;
my $groups   = undef;

GetOptions("kggfile=s"   =>   \$kggfile,
	   "datafile=s"  =>   \$datafile,
	   "groups=s"    =>   \$groups);

# get average exp per cluster

print "Calculating average profile per cluster ... ";
my $todo = "perl /Users/olivier/PERL_SCRIPTS/calculate_average_profile_per_group.pl $kggfile $datafile > $kggfile.avgexp";
system($todo);
print "Done.\n";

# use AOV to find clusters

print "Run ANOVA to find clusters that are differentially expressed ... ";

my @a_g   = split /\,/, $groups;
my $txt_g = "";
my $i     = 0;
my @a_s   = ();
my $sum   = 0;
foreach my $t (@a_g) {
  $txt_g = "rep($i,$t)";
  push @a_s, $txt_g;
  $sum += $t;
  $i++;
}

my $txt = "m <- read.csv(\"$kggfile.avgexp\", row.names=1, header=T, sep=\"\\t\")\n";
$txt .=   "groups <- c(" . join(", ", @a_s) . ")\n";
$txt .=   "pv <- as.matrix(apply(m, 1, function(x) { anova( aov( x ~ as.factor(groups) ) )\$'Pr(>F)'[1] } ))\n";
$txt .=   "write.table(rownames(m)[ pv < 0.01 / $sum ], file=\"$kggfile.avgexp.anova\", sep=\"\t\", row.names=F, quote=F)\n";

my $tmpfile = Sets::getTempFile("/tmp/Rscript");
open OUT, ">$tmpfile" or die "Cannot open $tmpfile for writing.\n";
print OUT $txt;
close OUT;
system("R CMD BATCH $tmpfile");

print "Done.\n";


print "Merge non-informative clusters ... ";
$todo = "perl ~/PERL_SCRIPTS/expression_merge_nonsignificant_clusters.pl $kggfile.avgexp.anova $kggfile > $kggfile.dereg";
system($todo);
print "Done.\n";

print "Create average expression profiles for these clusters only ... ";
$todo =  "perl /Users/olivier/PERL_SCRIPTS/calculate_average_profile_per_group.pl $kggfile.dereg $datafile > $kggfile.dereg.avgexp";
system($todo);
print "Done.\n";

print "Draw heatmap ... ";
$todo = "draw_expression_heatmap.pl --matrix=$kggfile.dereg.avgexp --draw=open --h=30 --clustrows=1";
system($todo);
print "Done\n";

