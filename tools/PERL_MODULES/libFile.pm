package      libFile;
require      Exporter;
@ISA       = qw(Exporter);
@EXPORT    = qw(append getcolumn get_hash get_exist_hash open_gz gzip_file 
		gunzip_file exist_zip file_access get_hash_hash);

# get hash of hashes...
sub get_hash_hash {
	my ($file, $col_key, $col1, $col2)= @_;

	if ( ! $file or ! -s $file ) {
	print STDERR "*** ERROR: input file $file not found\n";
	return undef;
	}

	$col_key= 0 if (! defined $col_key);
	$col1= 1 if (! defined $col1);
	$col2= 2 if (! defined $col2);

	open my $fin, $file or die $!;
	my ($line, @tmp);
	my %output= ();
	while ($line = <$fin>){
	chomp $line;
	if ($line =~ /^#|^\/\//){
	next;
	}
	@tmp= split(/\t/, $line);
	if (defined $tmp[$col_key] && defined $tmp[$col1] && defined $tmp[$col2]){ 
	# print "adding $tmp[$col_key] -> {$tmp[$col1]} -> $tmp[$col2]\n";
	$output{$tmp[$col_key]}{$tmp[$col1]}= $tmp[$col2];
	}
	}
	close $fin;
	return \%output;
}

sub append {
    my ($file1,$file2,$del_comment) = @_;
    my ($fh1,$fh2,$line) ;

    open ($fh1,">> $file1") or 
	return (0, "cannot open $file1 for appending:$!");
    open ($fh2,$file2) or 
	return (0, "cannot open $file2 for reading:$!");
    
    while ($line=<$fh2>) {
	next if ( $del_comment and $line =~ /^\#/ );
	print $fh1 $line;
    }
    close $fh2;
    close $fh1;
    return (1,"ok");
}

sub exist_zip {
    my ( $file ) = @_;
    my ( $file_zip, $file_unzip);

    $file_unzip = $file;
    if ( $file =~ /\.gz$/ ) {
	$file_unzip =~ s/\.gz$//;
	$file_zip = $file;
    } else {
	$file_zip = $file.'.gz';
    }
    if ( -s $file_zip or -s $file_unzip ) {
	return 1;
    } else {
	return 0;
    }
}

sub getcolumn {
    my ( $file, $column, $delim, $keep_comment ) = @_ ;
    $column = 0 if ( ! defined $column );
    $delim = "\t" if ( ! defined $delim );
 
    my ( $list,$fh,$line,@fields );

    return undef if ( ! -s $file or ! -r $file );

    $list = [ ];
    open ($fh,$file) or die "cannot open $file:$!";
    while ($line=<$fh>) {
	chomp $line;
	next if ( $line=~ /^\#/ and ! $keep_comment );
	next if ( $line !~ /\w+/ );

	@fields = split /$delim/, $line;
	push @$list, $fields[$column] if ( defined $fields[$column] );
    }
    close $fh;

    return $list;
}	


# get_hash:
# input: 
#	filename	duh.
#	column_key	specify column to read as key
#	column_value	specify column to read as value	
#
# commented by TS
sub get_hash {
    my ( $file,$column_key,$column_value ) = @_;
    my ( $fh,$line,@tmp, $key,$value );

    if ( ! $file or ! -s $file ) {
	print STDERR "*** ERROR: input file $file not found\n";
	return undef;
    }

    $column_key = 0 if ( ! defined $column_key );
    $column_value = 1 if ( ! defined $column_value );

    my $hash = { };


    open ($fh,$file) or die "cannot open $file:$!";
    while ( $line = <$fh> ) {
	next if ( $line !~ /\w/ );
	next if ( $line =~ /^\#/ );
	chomp $line;
	@tmp = split /\t/,$line;
	next if ( ! defined $tmp[$column_key] or ! defined $tmp[$column_value] );
	$key = $tmp[$column_key];
	$value = $tmp[$column_value];
	#do not remove spaces using lines below, it effects genome processing
	#$key=~s/\s+//g;
	#$value=~s/\s+//g;
	$hash->{$key} = $value;
    }
    close $fh;

    return $hash;
}

sub get_exist_hash {
    my ( $file,$column_key,$value ) = @_;
    my ( $fh,$line,@tmp, $key );

    if ( ! $file or ! -s $file ) {
	print STDERR "*** ERROR: input file not found\n";
	return undef;
    }

    $column_key = 0 if ( ! defined $column_key );
    $value = 1 if ( ! defined $value );

    my $hash = { };

    open ($fh,$file) or die "cannot open $file:$!";
    while ( $line = <$fh> ) {
	next if ( $line !~ /\w/ );
	next if ( $line =~ /^\#/ );
	chomp $line;
	@tmp = split /\t/,$line;
	next if ( ! defined $tmp[$column_key] );
	$key = $tmp[$column_key];
	$hash->{$key} = $value;
    }
    close $fh;

    return $hash;
}

sub open_gz {
    my ( $file ) = @_;
    my ( $file_zip, $file_unzip, $fileopen );

    undef $fileopen;

    $file_unzip = $file;
    if ( $file =~ /\.gz$/ ) {
	$file_unzip =~ s/\.gz$//;
	$file_zip = $file;
    } else {
	$file_zip = $file.'.gz';
    }
    if ( -s $file_unzip ) {
	$fileopen = "< $file_unzip ";
    } elsif ( -s $file_zip ) {
	$fileopen = " gunzip -c $file_zip | ";
    }

    return $fileopen;
}

sub gunzip_file {
    my ( $file, $force ) = @_;
    my ( $file_zip, $file_unzip);

    $file_zip = $file;
    if ( $file =~ /\.gz$/ ) {
	$file_unzip =~ s/\.gz$//;
    } else {
	$file_unzip = $file;
	$file_zip = $file.'.gz';
    }

    return 0 if ( ! -s $file_zip );
    if ( $force ) {
	system "gunzip -f $file_zip";
	return 2;
    }

    if ( -s $file_unzip ) {
	if ( -M $file_unzip > -M $file_zip ) {
	    system "gunzip -f $file_zip";
	} else {
	    unlink $file_zip;
	}
	return 2;
    } else {
	system "gunzip $file_zip";
    }

    return 1;
}

sub gzip_file {
    my ( $file, $force ) = @_;
    my ( $file_zip, $file_unzip);

    $file_unzip = $file;
    if ( $file =~ /\.gz$/ ) {
	$file_unzip =~ s/\.gz$//;
	$file_zip = $file;
    } else {
	$file_zip = $file.'.gz';
    }

    return 0 if ( ! -s $file_unzip );
    if ( $force ) {
	system "gzip -f $file_unzip";
	return 2;
    }

    if ( -s $file_zip ) {
	if ( -M $file_zip > -M $file_unzip ) {
	    system "gzip -f $file_unzip";
	} else {
	    unlink $file_unzip;
	}
	return 2;
    } else {
	system "gzip $file_unzip";
    }

    return 1;
}


sub file_access {
  my ($file,$days) =@_;#returns true if the file was accessed within the last $days
  if(!-s $file){
    print STDERR "input file $file not found\n";
    return;
  }
  my($ctime)=time;
  my($atime)=`stat --format=%Z $file`;#file access time
  chomp($atime);
  my($d)=($ctime-$atime)/(24*3600);
  if($d<$days){#was accessed in the last $days
    return 1;
  }
  else{
    return 0;
  }
}

1;
