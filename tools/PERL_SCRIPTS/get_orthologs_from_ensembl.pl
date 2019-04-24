use lib qw(/home/olly/PERL_MODULES);
use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my $i = 0;
my %H_UBRH = ();
my %H_MBRH = ();
my %H_RHS  = ();
my %H_BEST = ();

foreach my $r (@$a_ref) {
    if ($i == 0) {
	$i++; next;
    }
    if ($r->[2] eq "UBRH") {
	push @{ $H_UBRH{ $r->[0] } }, $r->[1];
	$H_BEST{ $r->[1] } = 1;
    }
    if ($r->[2] eq "MBRH") {
	push @{ $H_MBRH{ $r->[0] } }, $r->[1];
    }
    if ($r->[2] eq "RHS") {
	push @{ $H_RHS{ $r->[0] } }, $r->[1];
    }

    $i++;
}

my @s1 = keys(%H_UBRH);
my @s2 = keys(%H_MBRH);
my @s3 = keys(%H_RHS);


my $a_ref_genes_1 = Sets::getUnionSet(\@s1, \@s2);
my $a_ref_genes_2 = Sets::getUnionSet($a_ref_genes_1, \@s3);

my %H_TAKEN = ();

foreach my $k (@$a_ref_genes_2) {
    #print "$k => ". join("\t", @{ $H{$k} }); print "\n";

    if (defined($H_UBRH{$k})) {
	print "$k\t$H_UBRH{$k}->[0]\n";
    } elsif (defined($H_MBRH{$k})) {
	$H_MBRH{$k} = Sets::shuffle_array($H_MBRH{$k});
	foreach my $r (@{$H_MBRH{$k}}) {
	    if (!defined($H_TAKEN{ $r } )) {
		print "$k\t$r\n"; $H_TAKEN{ $r } = 1;
		last;
	    }
	}
	#print "$k\t$H_MBRH{$k}->[0]\n";
    }  elsif (defined($H_RHS{$k})) {
	$H_RHS{$k} = Sets::shuffle_array($H_RHS{$k});
	foreach my $r (@{$H_RHS{$k}}) {
	    if (!defined($H_BEST{ $r } ) && !defined($H_TAKEN{ $r } )) {
		print "$k\t$r\n"; $H_TAKEN{ $r } = 1;
		last;
	    }
	}
    } 
}
