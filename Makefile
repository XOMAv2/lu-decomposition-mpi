CC=mpicc
CFLAGS=-Wall -Wextra -pedantic -Werror

# $@ - цель
# $^ - зависимости цели
# $< - первый элемент из списка зависимостей

lu_decomposition: main.o matrix.o
	$(CC) -o $@ $^ $(CFLAGS)

main.o: main.c matrix.h
	$(CC) -c -o $@ $< $(CFLAGS)

%.o: %.c %.h
	$(CC) -c -o $@ $< $(CFLAGS)

.PHONY: clean

clean:
	rm -rf *.o lu_decomposition
