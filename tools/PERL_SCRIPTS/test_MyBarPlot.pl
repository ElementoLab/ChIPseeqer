BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use MyBarPlot;
use Table;

my $ta = Table->new;
$ta->loadFile($ARGV[0], 1, 1);
my $a_ref = $ta->getArray();
my $a_col = $ta->getColumn(0);

my $ba = MyBarPlot->new;
$ba->setMatrix($a_ref);
$ba->setRowMotif(1);
$ba->setRowNames($ta->getRowNames());
$ba->setHeader($ta->getHeader());

my $ta2 = Table->new;
$ta2->loadFile("$ARGV[0].pv", 1, 1);
my $a_ref_pv = $ta2->getArray();

$ba->setLabels($a_ref_pv);
$ba->setBarColor("red");
$ba->setXRight(200);
$ba->setMidLines($a_col);
$ba->setMinMax("0.0",0.5);
$ba->draw;
$ba->setOutputFileName("$ARGV[0].barplot.eps");

$ba->output();
$ba->pdfify();


