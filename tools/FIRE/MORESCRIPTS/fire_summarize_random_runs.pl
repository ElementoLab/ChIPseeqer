use Sets;

my $target_dir   = "$ARGV[0]\_RANDOM";

my $a_ref = Sets::getFiles("$target_dir/*.txt");

my $dna = 0;
my $rna = 0;
my $cnt_dna = 0;
my $cnt_rna = 0;
foreach my $f (@$a_ref) {
  my $fn = Sets::filename($f);
  my $ff = "$f\_FIRE";
  if (-e "$ff/DNA/$fn.summary") {
    $dna += Sets::nbLinesInFile("$ff/DNA/$fn.summary");
    $cnt_dna ++;
  }
  if (-e "$ff/RNA/$fn.summary") {
    $rna += Sets::nbLinesInFile("$ff/RNA/$fn.summary");
    $cnt_rna ++;
  }
}

my $cdna = sprintf("%4.3f", $dna / $cnt_dna);
my $crna = sprintf("%4.3f", $rna / $cnt_rna);

print "DNA: $cdna ($dna, $cnt_dna), RNA:$crna ($rna, $cnt_rna)\n"; 
