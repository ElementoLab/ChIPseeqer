#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";


use Table;
use Sets;
use strict;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my %IND = ();
my %H   = ();

foreach my $r (@$a_ref) {
    if ($r->[2] eq "CDS") {

	my ($ci) = $r->[8] =~ /Parent=(.+?)$/;
	
	#print "Pa=$ci\n";

	$IND { $ci } -> [0] = $ci;
	$IND { $ci } -> [1] = $r->[0];
	$IND { $ci } -> [2] = (defined($IND { $ci } -> [2])?Sets::min($IND { $ci } -> [2], $r->[3]):$r->[3]);
	$IND { $ci } -> [3] = Sets::max($IND { $ci } -> [3], $r->[4]);
	$IND { $ci } -> [4] = ($r->[6]eq'+'?1:-1);

    } elsif ($r->[2] eq "mRNA") {

	my ($ci) = $r->[8] =~ /ID=(.+?)\;/;	
	my ($al) = $r->[8] =~ /Alias=(.+?)$/;	

	$H{$ci}  = $al;

	#print "ID=$ci\n";

	$IND { $ci } -> [5] = (defined($IND { $ci } -> [5])?Sets::min($IND { $ci } -> [5], $r->[3]):$r->[3]);
	$IND { $ci } -> [6] = Sets::max($IND { $ci } -> [6], $r->[4]);
	
    }
}



foreach my $r (keys(%IND)) {
  $IND{$r}->[0] = $H{ $IND{$r}->[0] };
  print join("\t", @{ $IND{$r} }); print "\n";
}
