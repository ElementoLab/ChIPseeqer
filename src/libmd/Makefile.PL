# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.

use strict;
use ExtUtils::MakeMaker;
use Config;
use File::Copy;

my $arch = $Config{archname};
my $is_win32 = ($arch =~ /MSWin32/i);
my $is_387 = ($arch =~ /^i\d+-linux/i);
my $is_sol = ($arch =~ /sun4-solaris|sparc/i);
my $is_dar = ($arch =~ /darwin/i);
my $is_cyg = ($arch =~ /cygwin/i);
my $is_vms = ($arch =~ /vms/i);

unless ($is_387) {
    my $from = $is_win32 ? 'setpmsvc.c.win32' : 'setprec.c.unix';
    my $to = 'setprec.c';
    unless (-e $to) {
        copy($from, $to) or die "Cannot copy $from to $to: $!";
    }
}

unless ($is_387 or $is_sol) {
    my $from = 'sqrt.c.src';
    my $to = 'sqrt.c';
    unless (-e $to) {
        copy($from, $to) or die "Cannot copy $from to $to: $!";
    }
}

opendir(DIR, '.') or die "Cannot opendir '.': $!\n";
my @objs = map {s/\.c$/.o/; $_} grep { /\.c$/ } readdir DIR;
closedir DIR;
my %objs = map {$_ => 1} @objs;
foreach (qw(setprec.o sqrt.o)) {
    push @objs, $_ unless $objs{$_};
}

my %opts = (
            NAME      => 'Math::Cephes::libmd',
            VERSION_FROM => './Cephes.pm', # finds $VERSION
	    OBJECT => join(' ', @objs),
            SKIP      =>  [ qw( dynamic test ) ] ,
            LINKTYPE  =>  'static',
            clean => {FILES => 'libmd$(LIB_EXT)'},
	   );

#print "$Config{ccflags}\n";

$opts{CCFLAGS} = $Config{ccflags} . ' -Wall -fno-builtin ' 
    if $Config{gccversion};
$opts{CCFLAGS} = $Config{ccflags} . ' -arch x86_64 ' if ($Config{ccflags} =~ /arch ppc/);

WriteMakefile(%opts);

package MY;
sub c_o {
    my $self = shift;
    my $c_o = $self->SUPER::c_o(@_);
   if ($is_387) {
      $c_o .= << 'END';

sqrt.o: sqrtelf.387
	$(AS) -o sqrt.o sqrtelf.387

setprec.o: setprelf.387
	$(AS) -o setprec.o setprelf.387

END
   }
    elsif ($is_sol) {
       $c_o .= << 'END';

sqrt.o: sqrt.spa
	$(AS) -o sqrt.o sqrt.spa

END
   }
   else {

   }
   return $c_o;
}
  
sub post_constants {
    my $postconstant = <<'END';
INST_STATIC = libmd$(LIB_EXT)
END

  return $postconstant;
}

sub top_targets {
    my $top_targets = <<'END';
all    ::    static
	@$(NOOP)

static ::    libmd$(LIB_EXT)
	@$(NOOP)

config ::

test :	

END
   return $top_targets;
}
