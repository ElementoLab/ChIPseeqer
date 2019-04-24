BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Table;
use strict;

my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my %H  = ();
my %GO = ();
foreach my $r (@$a_ref) {
  my $n = shift @$r;
  foreach my $s (@$r) {
    #die "weird. H{$n}{$s} already seen\n" if (defined($H{$n}{$s}));
    $H{$n}{$s} = 1;
    $GO{$s}    ++;
  }  
}

$ta->loadFile($ARGV[1]);
my $h_ref_go = $ta->getIndexKV(0,1);

my @genes = keys(%H);
my @cats  = keys(%GO);


foreach my $c (@cats) {

  next if ($GO{$c} < 50);

  my $cc = $c;
  $cc =~ s/\:/\_/;
  my $gg = $h_ref_go->{$c};
  $gg =~ s/\ \(.+?\)//g;
  $gg =~ s/[\ \'\,\;\/\\]/\_/g; 
  $cc .= "\_$gg";
  print "\t$cc";
}
print "\n";

foreach my $g (@genes) {
  print "$g";
  foreach my $c (@cats) {
    
    next if ($GO{$c} < 50);

    if (defined($H{$g}{$c})) {
      print "\t1";
    } else {
      print "\t0";
    }
  }
  print "\n";
}
