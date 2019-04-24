my $infile = $ARGV[0];

my $todo   = undef;
if ($infile =~ /\.gz$/) {
  $todo = "gunzip $infile";
  system($todo) == 0 or die "Cannot exec $todo\n";
  $infile =~ s/\.gz$//;
}

$todo = "perl /panda_scratch_miro/ole2001/STAU2/table2fasta.pl $infile > $infile.fa";
system($todo) == 0 or die "Cannot exec $todo\n";
 
$todo = "bwa aln -t 4 /home/ole2001/PROGRAMS/ChIPseeqer-1.0/DATA/RepMask3.2.7_annotation_numberedIdx $infile.fa > $infile.fa.repeats";
system($todo) == 0 or die "Cannot exec $todo\n";

$todo = "bwa samse /home/ole2001/PROGRAMS/ChIPseeqer-1.0/DATA/RepMask3.2.7_annotation_numberedIdx $infile.fa.repeats $infile.fa > $infile.fa.repeats.sam";
system($todo) == 0 or die "Cannot exec $todo\n";

my $name = $infile;
$name =~ s/\_combined\_mastertable//;

$todo = "perl /home/ole2001/PROGRAMS/ChIPseeqer-1.0/SCRIPTS/CountRepeatMatches.pl --samfile=$infile.fa.repeats.sam --table=$infile --uniq=1 --name=$name > $infile.fa.repeats.sam.counts.uniq ";
system($todo) == 0 or die "Cannot exec $todo\n";

$todo = "perl /home/ole2001/PROGRAMS/ChIPseeqer-1.0/SCRIPTS/CountRepeatMatches.pl --samfile=$infile.fa.repeats.sam --table=$infile --uniq=1 --name=$name --t2c=1 > $infile.fa.repeats.sam.counts.uniq.t2c ";
system($todo) == 0 or die "Cannot exec $todo\n";

$todo = "perl /home/ole2001/PROGRAMS/ChIPseeqer-1.0/SCRIPTS/CountRepeatMatches.pl --samfile=$infile.fa.repeats.sam --table=$infile --uniq=0 --name=$name > $infile.fa.repeats.sam.counts ";
system($todo) == 0 or die "Cannot exec $todo\n";

$todo = "perl /home/ole2001/PROGRAMS/ChIPseeqer-1.0/SCRIPTS/CountRepeatMatches.pl --samfile=$infile.fa.repeats.sam --table=$infile --uniq=0 --name=$name --t2c=1 > $infile.fa.repeats.sam.counts.t2c ";
system($todo) == 0 or die "Cannot exec $todo\n";

