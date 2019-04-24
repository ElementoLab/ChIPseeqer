#
#  get a weblogo out of an AlignACE file
#

use lib qw(/home/olly/PERL_MODULES);
use MotifLibrary;
use Motif;
use Sets;

 



my $ml = MotifLibrary->new;
	    
$ml->loadAlignACEOutputFile($ARGV[1]);

my $mo = $ml->getOneMotif($ARGV[0]);


my $f = Sets::getTempFile("/tmp/gugu");
$mo->writeSites($f);

my $todo = "/home/olly/PERL_MODULES/PROGRAMS/weblogo/seqlogo -f $f -c -Y -n";
system($todo);
unlink $f;

