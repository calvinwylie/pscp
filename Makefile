CC=gcc
CFLAGS=-lglpk

.PHONY: clean

test: src/test.c
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm test