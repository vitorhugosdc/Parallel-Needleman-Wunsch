CC = mpicc
CFLAGS = 
TARGETS = sequencial threads mpi

all: $(TARGETS)

sequencial: sequencial.c
	$(CC) sequencial.c -o sequencial

threads: threads.c
	$(CC) $(CFLAGS) threads.c -o threads

mpi: mpi.c
	$(CC) $(CFLAGS) mpi.c -o mpi

clean:
	rm -f $(TARGETS)

.PHONY: all clean
