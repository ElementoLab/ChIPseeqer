#print 1 << 2;
for (my $i=0; $i<10; $i++) {
  if ($ARGV[0] & (1 << $i)) {
    print "1";
  } else {
    print "0";
  }

}
