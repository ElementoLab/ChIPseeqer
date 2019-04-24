open IN, $ARGV[0] or die "Cannot open $ARGV[0]\n";
while (my $l = <IN>) {
  chomp $l;
  my @a = split /\t/, $l, -1;

  my $desc = $a[8];

  $desc =~ s/gene_id/g_id/;
  $desc =~ s/transcript_id.+$//;
  
  print "$l$desc\n";

}
close IN;

