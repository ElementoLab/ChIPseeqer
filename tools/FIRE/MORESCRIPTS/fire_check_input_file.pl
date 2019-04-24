
use Getopt::Long;
use Table;
use Sets;
use strict;




my $expfile = undef;
my $exptype = undef;
GetOptions ('expfile=s'             => \$expfile,
	    'exptype=s'             => \$exptype);


my $ta = Table->new;
$ta->loadFile($ARGV[0]);
my $a_ref = $ta->getArray();

my $r = shift @$r; 
if ($r->[1] =~ /^[\d\.]*/) {
  print "WARNING: your file might not contain a header line.\n";
}

my %H = ();
my %V = ();
foreach my $r (@$a_ref) {
  if (defined($H{$r->[0]})) {
    die "Your files contains multiple rows with the same gene id. Please correct that before applying FIRE.\n";
  }
  $H{$r->[0]} ++;
  $V{$r->[1]} = 1;
} 


my @v = values ( %V );
my @v = sort { $a <=> $b } @v;
my $max = $v[$#v];
if (($exptype eq 'discrete') && (scalar(@v) != $max+1)) {
  die "Problem. Your discrete vector is missing some symbols.\n"
}
