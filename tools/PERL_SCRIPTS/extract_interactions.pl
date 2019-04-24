use lib qw(/home/elemento/PERL_MODULES);
use Table;


if (scalar(@ARGV) == 0) {
    die "Usage: .. interactions pattern type_inter\n";
}

my $ta = Table->new;
$ta->loadFile($ARGV[0]);

my $pat = $ARGV[1];
my $typ = $ARGV[2];

my $a_ref = $ta->getArray();

foreach my $r (@$a_ref) {

    next if (defined($typ) && ($r->[4] !~ /$typ/));
    
    #print join("\t", @$r); print "\n";

    my $s1 = $r->[2];
    my $s2 = $r->[3];

    my $pm = $r->[7];
    
    my @am = split /\;/, $pm, -1;
    shift @am; pop @am;

    $pm = join("/", @am);

    my $ty = $r->[4];


    my @a1 = split /\;/, $s1;
    my $n1 = undef;
    foreach my $t (@a1) {
	if ($t =~ /$pat/) {
	    $n1 = $t;
	}
    }
    my @a2 = split /\;/, $s2;
    my $n2 = undef;
    foreach my $t (@a2) {
	if ($t =~ /$pat/) {
	    $n2 = $t;
	}
    }
    $n1 =~ s/[a-z]$//;
    $n2 =~ s/[a-z]$//;
    
    if ($n1 && $n2) {
      push @{ $PMIDS{ $n1 }{ $n2 }{ $ty } }, @am; 
      #print "$n1\t$n2\t$pm\t$ty\n" 
    }
}

foreach my $n1 (keys( %PMIDS )) {
  
  foreach my $n2 (keys( %{ $PMIDS{ $n1 } } )) {
  
    foreach my $ty (keys( %{ $PMIDS{ $n1 }{ $n2 } } )) {
      
      my $pmids = join("/", @{ $PMIDS{ $n1 }{ $n2 }{ $ty } });

      print "$n1\t$n2\t$ty\t$pmids\n";
    
    }
    
  }
  
}
