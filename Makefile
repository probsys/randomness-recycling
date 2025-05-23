CFLAGS ?= -O3 -flto -march=native -Wno-unused-result

all: librr.a sample.out
	mkdir -p build/bin
	cp sample.out build/bin/sample_rr
	mkdir -p build/lib
	cp librr.a build/lib
	mkdir -p build/include
	cp *.h build/include
	$(MAKE) clean

%.o: %.c
	gcc $(CFLAGS) -c -o $@ $^

librr.a: types.o uniform.o binarysearch.o lookup.o alias.o aldr.o
	ar rcs $@ $^

%.out: %.c librr.a
	gcc $(CFLAGS) -o $@ $^ -lm

.PHONY: clean
clean:
	rm -rf *.a *.o *.out

test: all
	@echo "Running test..."
	./build/bin/sample_rr cdf 9000 1 1 2 3 2 | tr -d '\n' | tr ' ' '\n' | sort | uniq -c
	./build/bin/sample_rr lookup 9000 1 1 2 3 2 | tr -d '\n' | tr ' ' '\n' | sort | uniq -c
	./build/bin/sample_rr alias 9000 1 1 2 3 2 | tr -d '\n' | tr ' ' '\n' | sort | uniq -c
	./build/bin/sample_rr fldr 9000 1 1 2 3 2 | tr -d '\n' | tr ' ' '\n' | sort | uniq -c
	./build/bin/sample_rr aldr 9000 1 1 2 3 2 | tr -d '\n' | tr ' ' '\n' | sort | uniq -c
	cd examples && make
	./examples/example.out
