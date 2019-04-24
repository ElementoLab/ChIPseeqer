package      libBio;
require      Exporter;
@ISA       = qw(Exporter);
@EXPORT    = qw(load_fasta);


# make this a module
#	load_fasta(fname, column, nosuffix)
sub load_fasta {
	my $fname= shift;
	my $col= shift;
	my $nosuffix= shift;

	if (defined $nosuffix && $nosuffix =~ /^nos/){
		$nosuffix= 1;
	}else{
		$nosuffix= 0;
	}

	if (! defined $col){
		$col= 1;	# second column is id by default
	}

	my %seq= ();
	my $id;
	if (! -e $fname ){
		print STDERR "$fname does not exist!\n";
		return undef;
	}	
	open my $fh, $fname or die $!;
	while (my $line= <$fh>){
		chomp $line;
		if ($line =~ /^>/){
			$line= substr $line, 1;
			$id= (split(/\|| /, $line))[$col];
			if ($nosuffix){
				$id=~ s/\..+$//;
				# print "loading $id\n";
			}
		}else {	
			$line=~ s/ +//;
			$seq{$id}.= $line if ($id);
		}
	}
	close $fh;
	return \%seq;	
}
