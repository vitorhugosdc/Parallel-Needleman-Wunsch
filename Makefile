CC = cc
CFLAGS = -lpthread
TARGETS = sequencial paralela

all: $(TARGETS)

sequencial: sequencial.c
	$(CC) sequencial.c -o sequencial

paralela: paralela.c
	$(CC) $(CFLAGS) paralela.c -o paralela

clean:
	rm -f $(TARGETS)

.PHONY: all clean
