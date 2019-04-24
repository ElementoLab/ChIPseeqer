my $wget = 'wget';
if ($ENV{'LOGNAME'} eq 'tavlab') {
  $wget = '/Genomics/Users/elemento/usr/bin/wget';
}


my @species = ();
if ($ARGV[0] eq 'all') {
  @species = ('yeast', 'malaria', 'arabidopsis', 'drosophila', 'worm', 'mouse', 'human', 'pombe', 'ciona'); 
} else {
  @species = @ARGV;
}

  
foreach my $s (@species) {
  my $f = "$s\_data.zip";
  unlink $f if (-e $f); 
  print "Download $f\n";
  my $w = "http://tavazoielab.princeton.edu/FIRE/$f";  
  my $t = "$wget -O $f $w";
  print "Exec $t.\n";
  system($t);
  system("unzip -o $f");
}


