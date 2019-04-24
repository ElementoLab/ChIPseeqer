use Table;
use Sets;
use strict;

my $ta = Table->new;

$ta->loadFile($ARGV[0]);
my $a_ref1 = $ta->getArray();

$ta->loadFile($ARGV[1]);
my $a_ref2 = $ta->getArray();

my $txt1 = "";
my $cnt = 1;
foreach my $r (@$a_ref1) {

  my $wm = Sets::myre2wm($r->[0]);
  $txt1 .= "Motif $cnt\n$wm*********\n\n";
  
  $cnt ++;
}

my $txt2 = "";
my $cnt  = 1;
foreach my $r (@$a_ref2) {

  
  my $wm = Sets::myre2wm($r->[0]);
  
  if (defined($ARGV[2])) {
    #my $tmp = $r-;
    my $a   = Sets::shuffle_re($r->[0]);
    my $re  = join("", @$a);
    $wm     = Sets::myre2wm($re)  
    #print "$wm\n";
  }
  
  $txt2 .= "Motif $cnt\n$wm*********\n\n";
  
  $cnt ++;
}

open OUT1, ">tmp1";
print OUT1 $txt1;
close OUT1;

open OUT2, ">tmp2";
print OUT2 $txt2;
close OUT2;


my $todo = "alignace2004/CompareACE tmp1 tmp2 > tmp3";
system($todo);


$ta->loadFile("tmp3");
my $a_ref = $ta->getArray();

my %H1 = ();
my %V1 = ();
my %H2 = ();
my %V2 = ();
foreach my $r (@$a_ref) {

  #if ($r->[4] > 0.8) {
  #  print $a_ref1->[$r->[1]-1]->[0] . " matches " . $a_ref2->[$r->[3]-1]->[0] . "\n";
  #}

  if (!defined($V1{ $r->[1] }) || (defined($V1{$r->[1]}) && ($r->[4] > $V1{$r->[1]}))) {
    
    $V1{ $r->[1] } = $r->[4];
    $H1{ $r->[1] } = $r->[3];
  }

  if (!defined($V2{ $r->[3] }) || (defined($V2{$r->[3]}) && ($r->[4] > $V2{$r->[3]}))) {    
    $V2{ $r->[3] } = $r->[4];
    $H2{ $r->[3] } = $r->[1];
  }

}

foreach my $k1 (keys(%H1)) {
  
  my $k2 = $H1{$k1};

  next if ($V1{$k1} < 0.75);

  print $a_ref1->[$k1-1]->[0] . "\t" . $a_ref2->[$k2-1]->[0] . "\t" . $V1{$k1} . "\n";


  
}


unlink "tmp1";
unlink "tmp2";
unlink "tmp3";

