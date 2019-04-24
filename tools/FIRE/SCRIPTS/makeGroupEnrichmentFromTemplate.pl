

my $txt = undef;

my $usepm = 1;
if ($ARGV[0]) {  # if there is an argument, means use myhypergeom
  $usepm = 0;
  $txt = `perl -pi.bak -e 's/dummy2/$usepm/' GroupEnrichment.pm`;

} else {  # else use the 

  my $pwd = `pwd`; 
  $pwd    =~ s/[\r\n]//;
  my $out = `find $pwd/../modules -name \"Hypergeom.pm\"`;
  
  $out    =~ s/\/Hypergeom\.pm[\r\n]//;
  $out    =~ s/\//\\\//g;

  print "Path to Hypergeom.pm (determined by find): $out\n";
  my $txt = `perl -pi.bak -e 's/\\\/dummy1/$out/' GroupEnrichment.pm`;

}

