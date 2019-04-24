my $prog = shift @ARGV;

print "$prog : $prog.c dataio.o statistics.o information.o mi_library.o regexp.o sequences.o 
        \$(CC) \$(CFLAGS) -Wall -o $prog $prog.c dataio.o statistics.o information.o mi_library.o regexp.o sequences.o -lm
";
