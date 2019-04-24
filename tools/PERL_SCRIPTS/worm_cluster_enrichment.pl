require "/home/olly/CELEGANS/GO/Database.pm";


# what is there in this file ?
my @a = ();
while (my $l = <STDIN>) {

    chomp $l;

 
    next if ($l =~ /Cluster/);

    push @a, $l;

}

@a = map { $_ = "\"$_\"" } @a;

$s = join(", ", @a);


# pa utile now
#$sql = "SELECT distinct GO_ELEGANS.GO, GO_ELEGANS.NAME, TERM from GO_ELEGANS, GO_TERMS where GO_TERMS.GO = GO_ELEGANS.GO and GO_ELEGANS.NAME in ($s);";


$sql = "SELECT GO_ELEGANS.GO, count(distinct GO_ELEGANS.NAME) as COUNT, TERM from GO_ELEGANS, GO_TERMS where GO_TERMS.GO = GO_ELEGANS.GO and GO_ELEGANS.NAME in ($s) group by GO_ELEGANS.GO order by COUNT DESC;";

print "\n";


system("mysql -uroot ELEGANSHOUSE -e '$sql'");
