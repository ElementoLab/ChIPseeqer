BEGIN{ $home = `echo \$HOME`; chomp $home}

use lib "$home/PERL_MODULES";

use LWP::Simple;
use LWP;
use Sets;
use Log;

my $log = Log->new;

my $a_ref = Sets::readSet($ARGV[0]);


foreach my $r (@$a_ref) {
  
  my ($id) = $r =~ /CLSTR(\d{5})r1/;
  
  if (!$id) {
    $log->log("pb with $r\n");
  } else {
    my $txt = get("http://ghost.zool.kyoto-u.ac.jp/cgi-bin/txtget.cgi?$id");
    
    if ($txt !~ /Tailbud/) {
      $log->log("$id data incorrect ..\n");
    } else {
      open OUT, ">GHOST/$id.html";
      print OUT "$txt\n";
      close OUT;
    }
    sleep(1);
  }

  
}
