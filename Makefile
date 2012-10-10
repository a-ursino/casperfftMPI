# compiler Not Infiniband
CC = /opt/openmpi/bin/mpic++
# Compiler Infiniband Mellanoxl
ICC = /usr/mpi/gcc/openmpi-1.4.2/bin/mpic++
# looking for .h library path
HEADER = -I.
# Other compiler flags
CFLAGS = -Wall

OBJFILES := $(patsubst %.cpp,%.o, $(wildcard *.cpp))
HFILES := $(wildcard *.h)

3dfft:  $(OBJFILES)	
	$(CC) -o 3dfft.bin $(OBJFILES) $(CFLAGS)

%.o: %.cpp $(HFILES)
	$(CC) -c -o $@ $< $(CFLAGS)
clean:
	rm *.o 3dfft.bin

