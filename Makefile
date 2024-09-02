CC = cc
CFLAGS_THREADS = -lpthread
MPICC = mpicc
TARGETS = sequencial threads mpi

all: $(TARGETS)

sequencial: sequencial.c
	$(CC) sequencial.c -o sequencial

threads: threads.c
	$(CC) $(CFLAGS_THREADS) threads.c -o threads

mpi: mpi.c
	$(MPICC) mpi.c -o mpi

clean:
	rm -f $(TARGETS)

.PHONY: all clean
