package PBS;
use strict;

=head1 NAME

PBS

=head1 SYNOPSIS

   use PBS;

   my $pbs = PBS->new;

   # optional    
   $pbs->setPlatform($platform);
   $pbs->setQueue($queue);       
   $pbs->setMemory("4096Mb");

   # mandatory
   $pbs->setScriptName("script.pbs");
   $pbs->setWallTime("12:00:00");

   # now just add commands (that you would normally exexute manually)
   $pbs->addCmd("cd $pwd");
   $pbs->addCmd("export DYLD_LIBRARY_PATH=/Genomics/fafner/grid/users/elemento/usr/lib");
   $pbs->addCmd("echo \"DNA, remove duplicates, create $expfile_nodups_dna\"");

   # submit the script
   $pbs->submit; 

   # alternative: print the content of the script to be submitted
   $pbs->print;

   # alternative: just execute the script, without submitting (useful for debugging)
   $pbs->execute;


=cut


sub new {
    my ($class) = @_;

    my ($self) = {};

    $self->{SCRIPT_NAME} = 'script.pbs';
    $self->{WALLTIME}    = '24:00:00';
    $self->{MEMORY}      = undef; #'2000Mb';
    $self->{QUEUE}       = undef;
    $self->{ERRORS}      = undef;
    $self->{DEPJOBS}     = [];
    $self->{CMDS}        = [];
    $self->{PLATFORM}    = 'default';
    $self->{USEALLNODE}  = 0;
    bless $self;
    return $self;

}

sub getNumJobsForUser {
  my ($u) = @_;  
  my $txt = `qstat | grep $u`;
  my @a   = split /\n/, $txt;
  return scalar(@a);
}


sub useAllNode {
  my ($self, $f) = @_;
  $self->{USEALLNODE} = $f;
}

sub addDepJob {
  my ($self, $f) = @_;
  push @{ $self->{DEPJOBS} }, $f;
}


sub setQueue {
  my ($self, $f) = @_;
  $self->{QUEUE} = $f;
}


sub setWallTime {
  my ($self, $f) = @_;
  $self->{WALLTIME} = $f;
}


sub setPlatform {
  my ($self, $f) = @_;
  $self->{PLATFORM} = $f;
}


sub setMemory {
  my ($self, $f) = @_;
  $self->{MEMORY}   = $f;
}

sub setErrorFile {
  my ($self, $f) = @_;
  $self->{ERRORS}   = $f;
}

sub setScriptName {
  
  my ($self, $f) = @_;
  
  $self->{SCRIPT_NAME} = $f;

}


sub print {
  my ($self) = @_;
  
  my $txt = $self->_createText();
  
  print $txt;

}


sub addCmd {
  my ($self, $c) = @_;

  push @{  $self->{CMDS} }, $c;

  
}

sub _createText {
  my ($self) = @_;

  my $txt = "";
  $txt .= "#PBS -l mem=$self->{MEMORY}\n";
  $txt .= "#PBS -l walltime=$self->{WALLTIME}\n";


  $txt .= "#PBS -e $self->{SCRIPT_NAME}.e\n";
  $txt .= "#PBS -o $self->{SCRIPT_NAME}.o\n";
  
  #if (defined($self->{ERRORS})) {
  #}
  
  foreach my $c (@{ $self->{CMDS} }) {
    $txt .= "$c\n";
  }
    
  return $txt;

}

sub _writeScript {
  my ($self) = @_;
  
  my $txt = $self->_createText();
  open OUT, ">$self->{SCRIPT_NAME}";
  print OUT $txt;
  close OUT;

  

}

sub submitToGrid {
   my ($self) = @_;
   $self->_writeScript();
   system("qsub -cwd $self->{SCRIPT_NAME}");
}


sub submit {

  my ($self) = @_;

  $self->_writeScript();

  #if ($self->{PLATFORM} eq 'fafner') {
  system("chmod +x $self->{SCRIPT_NAME}");    
  #}

  my $todo = "qsub " ;
  if (($self->{PLATFORM} eq 'fafner') || ($self->{PLATFORM} eq 'tcluster'))  {
    $todo .= " -cwd ";
  }

  if ($self->{USEALLNODE} == 1) {
    $todo .= " -l nodes=1:ppn=2 ";
  }
  
  if (defined($self->{QUEUE})) {
    $todo .= " -q $self->{QUEUE} ";
  }
  
  # doc at  http://beige.ucs.indiana.edu/I590/node45.html
  #         http://www.arsc.edu/support/news/HPCnews/HPCnews320.shtml
  if (@{$self->{DEPJOBS}} > 0) {
    if ($self->{PLATFORM} eq 'tcluster') {
      $todo .= " -ac ";
    } else {
      $todo .= " -W ";
    }
    $todo .= " depend=afterany:" . join(":", @{$self->{DEPJOBS}});
  }
  
  $todo .= " $self->{SCRIPT_NAME}";

  print "$todo\n";

  my $out = `$todo`;
  $out =~ s/[\r\n]//g;

  if ($self->{PLATFORM} eq 'tcluster') {
    my ($realout) = $out =~ /Your\ job\ (\d+?)\ /;  
    $out = $realout;
  }

  return $out;
}


sub execute {
  my ($self) = @_;
   
  $self->_writeScript();

  system("sh $self->{SCRIPT_NAME}");
}

1;
