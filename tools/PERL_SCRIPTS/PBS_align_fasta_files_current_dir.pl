#!/usr/bin/perl
use lib "$ENV{HOME}/PERL_MODULES";

use Sets;
use PBS;
use strict;

my $a_ref = Sets::getFiles("*.fa");

foreach my $fa (@$a_ref) {

  while (PBS::getNumJobsForUser("ole2001") >= 50) {
    sleep(60);
  }
  my $fafa = $fa;
  $fafa =~ s/\.fa//;
  system("gr \"clustalw $fa\" $fafa");

}
