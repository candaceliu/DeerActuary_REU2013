
#define the suffixes
# .SUFFIXES:	.o .c .f 

#macro definitions
SHELL = /bin/sh
CC=mpic++
## CC=/opt/SUNWhpc/HPC8.2.1/gnu/bin/mpic++
CFLAGS=-O3 -Winline -finline-functions   -finline-limit=1200  -I/usr/local/include -Wall -std=c++0x
LIB=-lm -lpthread
LDEP=
DEP= 

.cpp.o:	
	$(CC) $(CFLAGS) $(CURLFLAGS) -c $<


all:	deerPopMilstein

deerPopMilstein: deerPopMilstein.o $(DEP) $(LDEP)
	$(CC) $(CFLAGS) -o $@ $@.o $(LIB)


clean:
	rm -f *.o *~ deerPopMilstein

