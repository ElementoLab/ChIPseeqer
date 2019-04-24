#!/usr/bin/perl
# takes as input new col, + STDIN ORF tab feature
use lib qw(/home/olly/PERL_MODULES);
use Database;;
use Getopt::Long;




GetOptions ('colname=s' => \$s_colname,
            'coltype=s' => \$s_coltype,
            'replace=s' => \$s_replace,
            'default=s' => \$s_default,
	    'after=s'   => \$s_after,
            'database=s' => \$s_database,
            'table=s'    => \$s_table,
            'index=s'    => \$s_index
            );

my $db = Database->new();
$db->connect($s_database);


if (!$s_colname || !$s_coltype) {    
    die "ins_col.pl --table=s --database=s --index=s --colname=s --coltype=s --replace=i --default=s --after=s\n";
}

if (!$s_colname && !$s_coltype) {
    die "please provide name and type";
}


if (!$s_replace) {

    my $s = "alter table $s_table add column $s_colname $s_coltype";

    if ($s_default) {

	if ($s_default eq 'NULL') {
	    $s .= " default NULL";
	} else {
	    
	    $s .= " default '$s_default'";
	    
	}
    }

    if ($s_after) {
	$s .= " after $s_after";
    }
    
    
    $db->query($s);
    
}


while (my $l = <STDIN>) {

    chomp $l;

    my ($p1, $p2) = split /\t/, $l;

    $p2 =~ s/\'/\ /g;

    my $s = "UPDATE $s_table set $s_colname = '$p2' where $s_index = '$p1'";
    
    #print "$s\n";

    $db->query($s);
    
}
