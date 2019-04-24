sub line_with_index {
    my $data_file   = shift;
    my $index_file  = shift;


    for ($line_number=0; $line_number<10000; $line_number+=100) {

      my $size;               # size of an index entry
      my $i_offset;           # offset into the index of the entry
      my $entry;              # index entry
      my $d_offset;           # offset into the data file
      
      $size = length(pack("L!", 0));
      $i_offset = $size * $line_number;
      seek($index_file, $i_offset, 0) or return;
      read($index_file, $entry, $size);
      $d_offset = unpack("L!", $entry);
      seek($data_file, $d_offset, 0);
      my $l = scalar(<$data_file>);
      
      my @a = split /\t/, $l, -1;

      print "$d_offset\t" . scalar(@a) . "\n";
    }
}

# usage:
$file = $ARGV[0];
open(FILE, "< $file")         or die "Can't open $file for reading: $!\n";
open(INDEX, "< $file.idx")
        or die "Can't open $file.idx for read/write: $!\n";
$line = line_with_index(*FILE, *INDEX, $ARGV[1]);

print $line;



