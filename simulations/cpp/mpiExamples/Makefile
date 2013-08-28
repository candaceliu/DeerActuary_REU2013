
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


all:	mpiBlockingCommunication mpiFileOperation mpiManagementExample mpiNonblockingCommunication

mpiBlockingCommunication: mpiBlockingCommunication.o $(DEP) $(LDEP)
	$(CC) $(CFLAGS) -o $@ $@.o $(LIB)

mpiFileOperation: mpiFileOperation.o $(DEP) $(LDEP)
	$(CC) $(CFLAGS) -o $@ $@.o $(LIB)

mpiManagementExample: mpiManagementExample.o $(DEP) $(LDEP)
	$(CC) $(CFLAGS) -o $@ $@.o $(LIB)

mpiNonblockingCommunication: mpiNonblockingCommunication.o $(DEP) $(LDEP)
	$(CC) $(CFLAGS) -o $@ $@.o $(LIB)


clean:
	rm -f *.o *~ mpiBlockingCommunication mpiFileOperation mpiManagementExample mpiNonblockingCommunication

