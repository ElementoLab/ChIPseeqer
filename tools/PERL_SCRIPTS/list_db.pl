#!/usr/bin/perl
use GDBM_File;

die "Please specify a file to list ..\n" 
    if !$ARGV[0];
tie(my(%db_nbfunc_text),'GDBM_File', $ARGV[0], &GDBM_READER, 0644);
reset(%db_nbfunc_text);
while (($k, $v) = each(%db_nbfunc_text)) {
    print "$k => $v\n";
}
