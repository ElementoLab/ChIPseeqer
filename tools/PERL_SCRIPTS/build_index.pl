my $file = $ARGV[0];

open(INDEX, "+>$file.idx") or die "Can't open $file.idx for read/write: $!\n";
open(NAMES, "+>$file.names") or die "blah blah ..\n";

open(FILE, "< $file") or die "Can't open $file for reading: $!\n";
build_index(*FILE, *INDEX);
close FILE;

open(FILE, "< $file") or die "Can't open $file for reading: $!\n";
build_names(*FILE, *NAMES);
close FILE;

close NAMES;
close INDEX;


sub build_index {
    my $data_file  = shift;
    my $index_file = shift;
    my $offset     = 0;
    while (<$data_file>) {
      print $index_file pack("L!", $offset);
      $offset = tell($data_file);      
    }
}


sub build_names {
    my $data_file  = shift;
    my $names_file = shift;

    while (my $l = <$data_file>) {
      my ($na) = $l =~ /^(.+?)\t/;
      print $names_file "$na\n";      
    }
}
