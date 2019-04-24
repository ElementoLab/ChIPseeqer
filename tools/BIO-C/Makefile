CFLAGS = -O3 -g -Wall
PROG = 
CC = cc
all : $(PROG)

$(PROG): $(PROG).c statistics.o information.o dataio.o hashtable.o interval_tree.o
	$(CC) $(CFLAGS)  -o $(PROG) $(PROG).c statistics.o information.o dataio.o hashtable.o interval_tree.o

mi_library.o : mi_library.c mi_library.h
	$(CC) $(CFLAGS)  -c mi_library.c 

statistics.o : statistics.c statistics.h
	$(CC) $(CFLAGS)  -c statistics.c 

sequences.o : sequences.c sequences.h
	$(CC) $(CFLAGS)  -c sequences.c 

dataio.o : dataio.c dataio.h
	$(CC) $(CFLAGS)  -c dataio.c 

information.o : information.c information.h
	$(CC) $(CFLAGS)  -c information.c 

regexp.o : regexp.c regexp.h
	$(CC) $(CFLAGS)  -c regexp.c 

prefix.o : prefix.c prefix.h
	$(CC) $(CFLAGS)  -c prefix.c 

hashtable.o : hashtable.c hashtable.h
	$(CC) $(CFLAGS)  -c hashtable.c

interval_tree.o : interval_tree.c interval_tree.h
	$(CC) $(CFLAGS)  -c hashtable.c

clean: 
	rm *.o

