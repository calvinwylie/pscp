CC=gcc-5
CFLAGS=-lglpk

SRC=src
BIN=bin

.PHONY: clean

test: $(BIN)/test.o $(BIN)/mt19937p.o
	mkdir -p $(BIN)/
	$(CC) $(CFLAGS) $^ -o $(BIN)/$@

portfolio: $(BIN)/portfolio.o $(BIN)/mt19937p.o
	mkdir -p $(BIN)/
	$(CC) $(CFLAGS) $^ -o $(BIN)/$@

$(BIN)/%.o: $(SRC)/%.c
	mkdir -p $(BIN)/
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm $(BIN)/*