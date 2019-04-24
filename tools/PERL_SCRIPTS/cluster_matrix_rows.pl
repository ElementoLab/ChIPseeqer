BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use strict;
use Matrix2EPS;
use AggloClust;
use Table;

my $h = 15;


my $ma = Matrix2EPS->new;
$ma->setOutputFileName("$ARGV[0].eps");
$ma->setMin(0);
$ma->setMax(5);
$ma->setMinText(" ");
$ma->setMaxText("Over-representation (log10)");
$ma->setVerbose(1);
$ma->setH($h);
$ma->setW(30);
$ma->setFont("Arial");
$ma->setAlgoClust("avg");

$ma->addColumnDesc($ARGV[1]);
$ma->setColMap("$ENV{HOME}/PROGRAMS/MYTREEVIEW/HEATMAPS/rb_cmap.csv");
$ma->drawRowNames("right");

$ma->setNeg(1);
$ma->loadMatrix($ARGV[0], 1, 1);

# $ma->setMatrix(\@M, \@motifs, \@regions);

$ma->clusterRows("euclidean");
$ma->clusterColumns("euclidean");


$ma->setXbase(50);
$ma->setYbase(250);
$ma->setXright(700);

#$ma->loadIclustPartition($ARGV[1]);
#$ma->reorderUsingPartition(1,1);

$ma->regexpColNames('s/\.txt\.ORF\.noamb$//');
$ma->regexpColNames('s/motif\_//');

$ma->colDrawMotifs(1);
$ma->setRownameFontSize($h);
$ma->setHeaderFontSize(30);

#$ma->loadColumnsPartition($ARGV[1]);


$ma->draw();
#$ma->addColumnClusters();
$ma->drawScale(20,380);

$ma->output();
$ma->pdfify();
