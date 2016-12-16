CC = gcc
FC = gfortran
CFLAGS = -O3
FFLAGS = -O3 -fimplicit-none

examples: basic_test basic_test_f77 indirect_writes indirect_writes_f77 tetrahedra_neighbours

basic_test: basic_test.c lplib3.c lplib3.h
	$(CC) $(CFLAGS) basic_test.c lplib3.c -lm -lpthread -o $@

tetrahedra_neighbours: tetrahedra_neighbours.c lplib3.c lplib3.h libmeshb7.c libmeshb7.h
	$(CC) $(CFLAGS) tetrahedra_neighbours.c lplib3.c libmeshb7.c -lm -lpthread -o $@

indirect_writes: indirect_writes.c lplib3.c lplib3.h libmeshb7.c libmeshb7.h
	$(CC) $(CFLAGS) indirect_writes.c lplib3.c libmeshb7.c -lm -lpthread -o $@

basic_test_f77: basic_test.f lplib3.c lplib3.h lplib3.ins
	$(CC) $(CFLAGS) -c lplib3.c
	$(FC) $(FFLAGS) -c basic_test.f
	$(FC) basic_test.o lplib3.o -lm -lpthread -o $@
	rm basic_test.o lplib3.o

indirect_writes_f77: indirect_writes.f lplib3.c lplib3.h lplib3.ins libmeshb7.c libmeshb7.h libmeshb7.ins
	$(CC) $(CFLAGS) -c lplib3.c libmeshb7.c -DF77API
	$(FC) $(FFLAGS) -c indirect_writes.f
	$(FC) indirect_writes.o lplib3.o libmeshb7.o  -lm -lpthread -o $@
	rm indirect_writes.o lplib3.o libmeshb7.o

clean:
	rm *.o basic_test basic_test_f77 indirect_writes indirect_writes_f77 tetrahedra_neighbours
