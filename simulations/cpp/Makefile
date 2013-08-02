
#define the suffixes
# .SUFFIXES:	.o .c .f 

#macro definitions
SHELL = /bin/sh
CC=g++
CFLAGS=-O3 -Winline -finline-functions   -finline-limit=1200  -I/usr/local/include
LIB=-lm
LDEP=
DEP= 

.cpp.o:	
	$(CC) $(CFLAGS) $(CURLFLAGS) -c $<


all:	deerPopMilstein

deerPopMilstein: deerPopMilstein.o $(DEP) $(LDEP)
	$(CC) $(CFLAGS) -o $@ $@.o $(LIB)


clean:
	rm -f *.o *~ deerPopMilstein

