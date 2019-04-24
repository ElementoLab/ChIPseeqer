
my $cntm = 0;

foreach my $f (@ARGV) {
  
  if ( -e $f) {

    open IN, $f; 
    my @a = <IN>;
    close IN;
    
    $cntm += scalar(@a);
  }

}

print "$cntm\n";
