# standard boilerplate 
CC	= cc
CFLAGS	+= -Wall -O3 -g -ggdb
LDFLAGS	= -lm

PROGKOH	= kohonen

.PHONY: all clean

all: $(PROGKOH)

$(PROGKOH):
	$(CC) $(CFLAGS) $(LDFLAGS) -o $(PROGKOH) $(PROGKOH).c
clean:
	rm -rf $(PROGKOH)
