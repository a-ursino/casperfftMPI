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

fft3d:  $(OBJFILES)	
	$(CC) -o fft3d.bin $(OBJFILES) $(CFLAGS)

%.o: %.cpp $(HFILES)
	$(CC) -c -o $@ $< $(CFLAGS)
clean:
	rm *.o fft3d.bin outputfft3d.o* outputfft3d.po*

