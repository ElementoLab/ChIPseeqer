CFLAGS += -O3 -g 

CC = gcc
all : page

page : page.c statistics.o dataio.o information.o mi_library.o hashtable.o
	$(CC) $(CFLAGS) -Wall -o page page.c statistics.o dataio.o information.o mi_library.o hashtable.o -lm $(LFLAGS)

mi_go_motif_calculator : mi_go_motif_calculator.c statistics.o dataio.o information.o mi_library.o
	$(CC) $(CFLAGS) -Wall -o mi_go_motif_calculator mi_go_motif_calculator.c statistics.o dataio.o information.o mi_library.o -lm $(LFLAGS) -O3

mi_analyze_matrix : mi_analyze_matrix.c statistics.o dataio.o information.o mi_library.o
	$(CC) $(CFLAGS) -Wall -o mi_analyze_matrix mi_analyze_matrix.c  statistics.o dataio.o information.o mi_library.o $(LFLAGS)

mi_library.o : mi_library.c mi_library.h
	$(CC) $(CFLAGS) -Wall -c mi_library.c 

statistics.o : statistics.c statistics.h
	$(CC) $(CFLAGS) -Wall -c statistics.c 

dataio.o : dataio.c dataio.h
	$(CC) $(CFLAGS) -Wall -c dataio.c 

information.o : information.c information.h
	$(CC) $(CFLAGS) -Wall -c information.c 

hashtable.o : hashtable.c hashtable.h
	$(CC) $(CFLAGS) -Wall -c hashtable.c 

clean: 
	rm -rf *.o *.dSYM page mi_go_motif_calculator mi_analyze_matrix
