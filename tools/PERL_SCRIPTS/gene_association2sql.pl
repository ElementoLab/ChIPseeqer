use lib qw(/home/olly/PERL_MODULES);
use GO;
use GO_func;
use Database;
use strict;
use DataFiles;
use strict;

my $df = DataFiles->new;

my $verbose = 1;

my $type = (defined($ARGV[1])?$ARGV[1]:"P");

my $go = GO->new;

if ($type eq "F") {
    #$go->loadOntology($df->get("GO_FUNCTION_ONTOLOGY"));
    $go->loadOntology("./function.ontology");
} elsif ($type eq "P") {
    #$go->loadOntology($df->get("GO_PROCESS_ONTOLOGY"));
    $go->loadOntology("./process.ontology");
} else {
    $go->loadOntology("./component.ontology");
    #$go->loadOntology($df->get("GO_COMPONENT_ONTOLOGY"));
}

open IN, $ARGV[0] or die "Cannot open tab file\n";

my %annot_go = ();

my $sp = $ARGV[2];

#
#  get all the GO categories associated with a single gene
#
while (my $l = <IN>) {
    chomp $l;

    my @a = split /\t/, $l;

    $a[4] =~ s/GO\://; $a[4] = int($a[4]);

    $a[10] =~ s/[\'\"]//g;
    $a[2]  =~ s/[\'\"]//g;
    my %h = (
	UID       => uc($a[1]),
	US        => $a[2],
	GOID      => $a[4],
	EVIDENCE  => $a[6],
	ASPECT    => $a[8],
	SYNONYM   => $a[10]
	);
    
    if ($sp eq "wb") {
	$h{UID} = ($h{SYNONYM} ne ""?$h{SYNONYM}:$h{US});
    }

    #next if ($h{US} ne 'SOC2_HUMAN');

    next if ($a[8] ne $type);

	
    next  if ($a[1] eq "");


    $a[1] = uc($a[1]);

    push @{ $annot_go{ $a[1] } }, \%h;
	    
	    
}


#
#  for each gene, get all the OID associated, and retrieve the parents
#
foreach my $uid (keys(%annot_go)) {

    #print ">>$uid\n" if ($verbose);
    
    # get all the GOID associated
    my @goids = ();
    foreach my $h ( @{$annot_go{ $uid }} ) {
	push @goids, $h->{GOID};
    }

    #print join(" ", @goids);
    #print "\n";
    
    # for each uid
    foreach my $h ( @{$annot_go{ $uid }} )  {

	my $thegoid = $h->{GOID};

	my @a_parents = @{$go->getParentVector( [ $h->{GOID} ] )};

	foreach my $p (@a_parents) {
	    
	    # process tops
	    next if ( ($p == 4) || ( $p == 3673 ) || ( $p == 8150) );
	    
	    # function tops
	    next if ( $p == 3674 );

	    # component tops
	    next if ( ( $p == 5575 ) || ( $p == 5623) );

	    if ($thegoid == $p) {
		$h->{ORIGINAL} = 1;
	    } else {
		$h->{ORIGINAL} = 0;
	    }

	    $h->{GOID} = $p;
	    
	    print Database::hash2insert($h, "GO_FULL_ANNOTATION");
	    print "\n";
	    
	}
    }
}


