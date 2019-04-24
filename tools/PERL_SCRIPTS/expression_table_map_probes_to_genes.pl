BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";
use Table;

my $ta = Table->new;
$ta->loadFile((defined($ARGV[1])?$ARGV[1]:"/home/elemento/PROJECTS/MBT_TRANSITION/DATA/probe_annotation_newcg.txt"));
my $h_ref = $ta->getIndexKV((defined($ARGV[2])?$ARGV[2]:0),(defined($ARGV[3])?$ARGV[3]:3));

use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my $cnt = 0;
foreach my $r (@$a_ref) {
  if ($cnt == 0) {
    print join("\t", @$r); print "\n";
  } else { 
  $r->[0] = $h_ref->{ $r->[0] } unless (($cnt == 0) && ($r->[0] !~ /\_at/));   

  if (defined($r->[0]) && ($r->[0] ne "NA")) {

    if ($r->[0] =~ /\//) {
      my @a = split /\//, $r->[0];
      foreach my $s (@a) {
	$r->[0] = $s;
	print join("\t", @$r); print "\n";
      }
    } else {
      print join("\t", @$r); print "\n";
    }
  }
 }
  $cnt ++;
}

