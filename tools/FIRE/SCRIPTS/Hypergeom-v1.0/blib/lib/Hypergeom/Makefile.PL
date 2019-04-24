use ExtUtils::MakeMaker;
$Verbose = 1;
WriteMakefile(
	NAME   => 'Hypergeom::hypergeom',
	SKIP   => [qw(all static static_lib dynamic dynamic_lib)],
	clean  => {'FILES' => 'hypergeom$(LIB_EXT)'},
);

sub MY::top_targets {
    '
all :: static

pure_all :: static

static ::       hypergeom$(LIB_EXT)

hypergeom$(LIB_EXT): $(O_FILES)
	$(AR) cr hypergeom$(LIB_EXT) $(O_FILES)
	$(RANLIB) hypergeom$(LIB_EXT)

';
}  
