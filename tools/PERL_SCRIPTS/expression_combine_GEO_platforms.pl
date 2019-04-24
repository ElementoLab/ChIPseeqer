my $gse = undef;
foreach my $f (@ARGV) {
  my @a = split /[\-\_]/, $f;
  system("getgeoplatform $a[1]");
  $gse = $a[0];
  my $txt = `perl ~/PERL_MODULES/SCRIPTS/show_matching_matrix_entries.pl $a[1].txt T01E8.5`;
  my @b   = split /[\n\r]+/, $txt;
  my %c   = ();
  foreach my $r (@b) {
    my @d = split /\t/, $r;
    $c{$d[1]} = 1; 
  }
  if (scalar(keys(%c)) != 1) {
    die "Problem with column of T01E8.5\n";
  } 

  my @e = keys(%c);
  
  system("sh ~/PROGRAMS/MIMOTIFS/TEMPORARY/GEO_process_series.sh $f $a[1].txt 0 $e[0]");

}


if (!defined($gse)) {
  die "gse not defined, cannot continue\n";
}
system("combine_geo_platforms \"*.rowavg\" $gse-combined.txt");

print "Cluster into :\n";
system("getnbclusters $gse-combined.txt")

