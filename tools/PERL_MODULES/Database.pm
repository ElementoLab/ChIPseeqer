package Database;
use strict;
use DBI;



sub new {    
    my $self  = {};
    $self->{DBH} = undef;
    $self->{VERBOSE} = 0;
    $self->{USER}    = "root";
    $self->{PASS}    = "";
    $self->{HOST}    = "localhost";
    
    bless($self);           # but see below
    return $self;    

}

#
#  
#
sub setVerbose {
    my ($self, $d) = @_;
    $self->{VERBOSE} = $d;
}

sub connect {
    
    my ($self, $db) = @_;

    $self->{DBH} = DBI->connect(    "DBI:mysql:database=$db;host=" . $self->{HOST}, 
				    $self->{USER}, 
				    $self->{PASS},
				    {'RaiseError' => 1}
				    );
    return $self->{DBH};
}


#
#  set the ids for the connection
#
sub setID {
    my ($self, $u, $p, $h) = @_;
    
    $self->{USER} = $u;
    $self->{PASS} = $p;
    $self->{HOST} = $h;
    
}

#
#  disconnect from the database
#
sub disconnect {
    my ($self) = @_;

     $self->{DBH}->disconnect; 
}


#
# recupere un tableau a partir d'une requete SQL
#
sub queryAllRecords {
        my ($self, $lstr_SQLquery) = @_;


	print "Executing: $lstr_SQLquery\n" if ($self->{VERBOSE} == 1);
        my $sth = $self->{DBH}->prepare($lstr_SQLquery);
	my $res = $sth->execute() or die("pb avec $lstr_SQLquery\n");
        my @rows = ();
        while (my $row = $sth->fetchrow_hashref) {
	    #print "pushing $row";
	    push(@rows, $row);
        }
        return @rows;
}


sub queryAllRecordsRef {
        my ($self, $lstr_SQLquery) = @_;


	print "Executing: $lstr_SQLquery\n" if ($self->{VERBOSE} == 1);
        my $sth = $self->{DBH}->prepare($lstr_SQLquery);
	my $res = $sth->execute() or die("pb avec $lstr_SQLquery\n");
        my @rows = ();
        while (my $row = $sth->fetchrow_hashref) {
	    #print "pushing $row";
	    push(@rows, $row);
        }
        return \@rows;
}

sub queryAllRecordsArrayRef {
        my ($self, $lstr_SQLquery) = @_;


	print "Executing: $lstr_SQLquery\n" if ($self->{VERBOSE} == 1);
        my $sth = $self->{DBH}->prepare($lstr_SQLquery);
	my $res = $sth->execute() or die("pb avec $lstr_SQLquery\n");
        my @rows = ();
        while (my $row = $sth->fetch) {
	    my @a_tmp = @$row;
	    # print "PUSH ($row)" . join("\t", @$row); print "\n";
	    push @rows, \@a_tmp;
        }


        return \@rows;
}




sub queryOneRecord {
        my ($self, $lstr_SQLquery) = @_;

	print "Executing: $lstr_SQLquery\n" if ($self->{VERBOSE} == 1);

        my $sth = $self->{DBH}->prepare($lstr_SQLquery);

        my $res = $sth->execute() or die("pb avec $lstr_SQLquery\n");
        my $row = $sth->fetchrow_hashref;

        return $row;

}

sub queryOneRecordArrayRef {
        my ($self, $lstr_SQLquery) = @_;

	print "Executing: $lstr_SQLquery\n" if ($self->{VERBOSE} == 1);

        my $sth = $self->{DBH}->prepare($lstr_SQLquery);

        my $res = $sth->execute() or die("pb avec $lstr_SQLquery\n");
        my $row = $sth->fetchrow_arrayref;

        return $row;

}


sub query 
{
        my ($self, $query) = @_;

	print "Executing: $query\n" if ($self->{VERBOSE} == 1);
	

        my $sth = $self->{DBH}->do($query) or die("pb avec $query\n");
        
}

#
#  return the last insert id
#
sub lastInsertId {
    my ($self) = @_;
    return $self->{DBH}->{mysql_insertid};
    
}


#
# transform a hash into a SQL insert
#
sub hash2insert {
    my ($h_ref, $table) = @_;
    
    my @p1 = ();
    my @p2 = ();

    reset(%$h_ref);
    while (my ($k, $v) = each(%$h_ref)) {
	#print "$k => $v\n";
	push @p1, $k;
	push @p2, $v;
    }

    map {$_ = "'$_'"} @p2;
    
    my $tmp = "insert into $table (" . 
	join(", ", @p1) . ") values (" . join(", ", @p2) . ");";

    return $tmp;
}


1;
