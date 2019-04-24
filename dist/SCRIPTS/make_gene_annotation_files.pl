#!/usr/bin/perl

use lib "$ENV{CHIPSEEQERDIR}";

use Getopt::Long;

# variables to store the arguments values
my $annotation		= undef;

# handling lack of arguments
if (@ARGV == 0) {
	die "Usage: make_gene_annotation_files.pl --annotation=FILE \n";
}

# processing command line options
GetOptions("annotation=s" => \$annotation );


#
# Make NM2ORF file
#
my $todo = "awk -F \"\t\" '{print \$2 \"\t\" \$13}' $annotation > $annotation.NM2ORF";
system($todo) == 0 or die "Cannot exec $todo\n";
if (-e "$annotation.NM2ORF") {
  print "Created $annotation.NM2ORF\n";
} else {
  print "Could not create\n";
}
#
# Make small annotation file
#
$todo = "awk -F \"\t\" '{print \$2 \"\t\" \$3 \"\t\" \$7 \"\t\" \$8 \"\t\" \$4 \"\t\" \$5 \"\t\" \$6}' $annotation > $annotation.new";
system($todo) == 0 or die "Cannot exec $todo\n";

$todo= "perl -pi -e 's/\\+/1/g' $annotation.new";
system($todo) == 0 or die "Cannot exec $todo\n";
$todo= "perl -pi -e 's/\-/-1/g' $annotation.new";
system($todo) == 0 or die "Cannot exec $todo\n";
if (-e "$annotation.new") {
  print "Created $annotation.new\n";
} else {
  print "Could not create\n";
}

#
# Create exons file
#
$todo= "perl $ENV{CHIPSEEQERDIR}/SCRIPTS/refGene2exons.pl $annotation > $annotation.EXONS.tmp";
system($todo) == 0 or die "Cannot exec $todo\n";
$todo= "awk -F\"\t\" '{print \"EX-\"\$1 \"\t\" \$2 \"\t\" \$3 \"\t\" \$4}' $annotation.EXONS.tmp > $annotation.EXONS";
system($todo) == 0 or die "Cannot exec $todo\n";
if (-e "$annotation.EXONS") {
  print "Created $annotation.EXONS\n";
} else {
  print "Could not create\n";
}
unlink("$annotation.EXONS.tmp");

#
# Create introns file
#
$todo= "perl $ENV{CHIPSEEQERDIR}/SCRIPTS/refGene2introns.pl $annotation > $annotation.INTRONS.tmp";
system($todo) == 0 or die "Cannot exec $todo\n";
$todo= "awk -F\"\t\" '{print \"I-\"\$0}' $annotation.INTRONS.tmp > $annotation.INTRONS";
system($todo) == 0 or die "Cannot exec $todo\n";
if (-e "$annotation.INTRONS") {
  print "Created $annotation.INTRONS\n";
} else {
  print "Could not create\n";
}
unlink("$annotation.INTRONS.tmp");

#
# Create GENEPARTS file
#
$todo= "cat $annotation.INTRONS >> $annotation.GENEPARTS";
system($todo) == 0 or die "Cannot exec $todo\n";
$todo= "cat $annotation.EXONS >> $annotation.GENEPARTS";
system($todo) == 0 or die "Cannot exec $todo\n";
if (-e "$annotation.GENEPARTS") {
  print "Created $annotation.GENEPARTS\n";
} else {
  print "Could not create\n";
}
#
# Create TSS_TES file
#
$todo= "awk -F\"\t\" '{print \$1 \"\t\" \$2 \"\t\" \$6 \"\t\" \$7 \"\t\" \$5}' $annotation.new > $annotation.TSS_TES";
system($todo) == 0 or die "Cannot exec $todo\n";
if (-e "$annotation.TSS_TES") {
  print "Created $annotation.TSS_TES\n";
} else {
  print "Could not create\n";
}
#
# Create oneperTSS file (used in ChIpseeqer2DensityMatrix)
#
$todo= "perl $ENV{CHIPSEEQERDIR}/SCRIPTS/refGene_retain_one_RefSeq_transcript_per_promoter.pl $annotation 1 > $annotation.oneperTSS";
system($todo) == 0 or die "Cannot exec $todo\n";
if (-e "$annotation.oneperTSS") {
  print "Created $annotation.oneperTSS\n";
} else {
  print "Could not create\n";
}
