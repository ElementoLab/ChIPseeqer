use lib "$ENV{HOME}/PERL_MODULES";
use PBS;

while (1) {
  my $p = PBS::getNumJobsForUser("ole2001");
  if ($ARGV[0] ne "") {
    print "$p jobs left.\n";
  }
  if ($p == 0) {
    last;
  }
  sleep(60);
}

system("echo \"nothing\" | mail -s \"Jobs Finished\" ole2001\@med.cornell.edu");
