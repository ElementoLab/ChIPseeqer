CFLAGS += -g -O3
CC = cc

PROGCOOC = EvaluateMotifColocalization

all : MyScanACE $(PROGCOOC)

MyScanACE : MyScanACE.c statistics.o dataio.o information.o sequences.o hashtable.o
	$(CC) $(CFLAGS) -Wall -o MyScanACE MyScanACE.c statistics.o dataio.o information.o sequences.o hashtable.o -lm

$(PROGCOOC) : $(PROGCOOC).c statistics.o dataio.o information.o sequences.o hashtable.o
	$(CC) $(CFLAGS) -Wall -o $(PROGCOOC) $(PROGCOOC).c statistics.o dataio.o information.o sequences.o hashtable.o -lm

mi_library.o : mi_library.c mi_library.h
	$(CC) $(CFLAGS) -Wall -c mi_library.c 

statistics.o : statistics.c statistics.h
	$(CC) $(CFLAGS) -Wall -c statistics.c 

hashtable.o : hashtable.c hashtable.h
	$(CC) $(CFLAGS) -Wall -c hashtable.c 

dataio.o : dataio.c dataio.h
	$(CC) $(CFLAGS) -Wall -c dataio.c 

information.o : information.c information.h
	$(CC) $(CFLAGS) -Wall -c information.c 

regexp.o : regexp.c regexp.h
	$(CC) $(CFLAGS) -Wall -c regexp.c 

prefix.o : prefix.c prefix.h
	$(CC) $(CFLAGS) -Wall -c prefix.c 

sequences.o : sequences.c sequences.h
	$(CC) $(CFLAGS) -std=c99 -Wall -c sequences.c 

clean: 
	rm -rf *.o *.dSYM MyScanACE $(PROGCOOC)

