# read lanes.txt
open IN, $ARGV[0] or die "Cannot open $ARGV[0]\n";
my @files = ();
my @labels = ();
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;
  my $prefpkm = undef;
  if ($a[0] =~ /^\d+$/) {    
    $prefpkm = "s_$a[0]_sequence.txt_tophat/transcripts.expr";
  } else {
    $prefpkm = "$a[0]/transcripts.expr";
  }
  
  my $todo = "columns.pl 0 5 < $prefpkm > $prefpkm.fpkm";
  system($todo) == 0 or die "Cannot $todo\n";
 
  push @files, "$prefpkm.fpkm";
  push @labels, $a[1];
}
close IN;

my $todo = "expression_concatenate_matrices.pl " . join(" ", @files);

system($todo) == 0 or die "Cannot $todo\n";

