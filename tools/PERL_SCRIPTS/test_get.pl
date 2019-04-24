use LWP::Simple;
use LWP;

$t = get("http://www.yeastgenome.org");

print $t;

