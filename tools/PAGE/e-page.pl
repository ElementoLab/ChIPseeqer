#!/usr/bin/perl
BEGIN{ $home = `echo \$HOME`; chomp $home}
use lib "$home/PERL_MODULES";

use Getopt::Long;
use strict;

my $expfile        = undef;
my $exptype        = undef;
my $profiles       = undef;
my $profiles_data  = undef;
my $profiles_names = undef;
my $maxnumprofiles = 10;
my $quantized      = 1;
my $epagedir       = "$ENV{HOME}/PROGRAMS/CONNECTIVITY_MAP";
my $suffix         = undef;
my $max_p          = 0.005;
my $verbose        = 1;
my $ebins          = undef;
my $outprofiles    = undef;

GetOptions("expfile=s"  =>  \$expfile,
	   "suffix=s"   =>  \$suffix,
	   "ebins=s"    =>  \$ebins,
	   "profiles=s" =>  \$profiles,
	   "outprofiles=s" =>  \$outprofiles,	   
	   "exptype=s"  =>  \$exptype);

if (defined($suffix)) {

  if ($suffix =~ /^pathway/) {
    $suffix = $profiles;
  }

  my $newexpfile = "$expfile.$suffix";
  system("cp $expfile $newexpfile");
  $expfile = $newexpfile;
}

#
# Making workspace directory
#
if (! -d "$expfile\_PAGE") 
{
    mkdir("$expfile\_PAGE") or die "couldn't make the directory: $?";
}


if ($profiles eq "cmap") {
  $profiles_data  = "$epagedir/DATA/ratioMatrix.txt.rowavg";
  $profiles_names = "$epagedir/DATA/cmap_instances_02.txt";
}


if ($exptype eq "continuous") {
  $quantized = 0;
}


my $todo = "$epagedir/e-page -expfile $expfile -datafile $profiles_data ";
$todo .= " -expnames $profiles_names -quantized $quantized ";
$todo .= " -maxnumprofiles $maxnumprofiles ";
$todo .= " -outprofiles $expfile\_PAGE/signif.profiles.txt -max_p $max_p ";
$todo .= " -outdata $expfile\_PAGE/signif.data.txt ";
$todo .= " -outfile $expfile\_PAGE/signif.out.txt ";

if (defined($ebins)) {
  $todo .= " -ebins $ebins ";
}

if ($verbose == 1) {
  print "$todo\n";
  $todo .= " -verbose 1 ";
}

system($todo);
