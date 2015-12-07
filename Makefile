CC=icc
CFLAGS=-fopenmp -I/share/apps/utils/include -L/share/apps/utils/lib -lglpk \
       -O3 -no-prec-div -xcore-avx2 -ipo

SRC=src
BIN=bin

.PHONY: clean

test: $(BIN)/test.o $(BIN)/mt19937p.o
	mkdir -p $(BIN)/
	$(CC) $(CFLAGS) $^ -o $(BIN)/$@

portfolio: $(BIN)/portfolio.o $(BIN)/mt19937p.o $(BIN)/rnglib.o $(BIN)/ranlib.o
	mkdir -p $(BIN)/
	$(CC) $(CFLAGS) $^ -o $(BIN)/$@

$(BIN)/%.o: $(SRC)/%.c
	mkdir -p $(BIN)/
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm $(BIN)/*