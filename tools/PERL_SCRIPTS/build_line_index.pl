my $file = $ARGV[0];
open(FILE, "< $file") or die "Can't open $file for reading: $!\n";
open(INDEX, "+>$file.idx") or die "Can't open $file.idx for read/write: $!\n";
build_index(*FILE, *INDEX);

sub build_index {
    my $data_file  = shift;
    my $index_file = shift;
    my $offset     = 0;

    while (<$data_file>) {
        print $index_file pack("N", $offset);
        $offset = tell($data_file);
    }
}
