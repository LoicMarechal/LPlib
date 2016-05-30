CC = /usr/local/bin/gcc
FC = /usr/local/bin/gfortran
CFLAGS = -Ofast
FFLAGS = -Ofast -fimplicit-none

examples: basic_test basic_test_f77 indirect_writes indirect_writes_f77 tetrahedra_neighbours #build_edges

basic_test: basic_test.c lplib3.c lplib3.h
	$(CC) $(CFLAGS) basic_test.c lplib3.c -lm -lpthread -o $@

basic_test_f77: basic_test.f lplib3.c lplib3.h lplib3.ins
	$(CC) $(CFLAGS) -c lplib3.c
	$(FC) $(FFLAGS) -c basic_test.f
	$(FC) basic_test.o lplib3.o -lm -lpthread -o $@
	rm basic_test.o lplib3.o

indirect_writes: indirect_writes.c lplib3.c lplib3.h libmesh6.c libmesh6.h
	$(CC) $(CFLAGS) indirect_writes.c lplib3.c libmesh6.c -lm -lpthread -o $@

indirect_writes_f77: indirect_writes.f lplib3.c lplib3.h lplib3.ins libmesh6.c libmesh6.h libmesh6.ins
	$(CC) $(CFLAGS) -c lplib3.c libmesh6.c 
	$(FC) $(FFLAGS) -c indirect_writes.f
	$(FC) indirect_writes.o lplib3.o libmesh6.o  -lm -lpthread -o $@
	rm indirect_writes.o lplib3.o libmesh6.o

tetrahedra_neighbours: tetrahedra_neighbours.c lplib3.c lplib3.h libmesh6.c libmesh6.h
	$(CC) $(CFLAGS) tetrahedra_neighbours.c lplib3.c libmesh6.c -lm -lpthread -o $@

build_edges: build_edges.c lplib3.c lplib3.h libmesh6.c libmesh6.h
ifdef NUMA
	$(CC) -O3 -funroll-loops -I . build_edges.c lplib3.c libmesh6.c -lm -lpthread -L . -lnuma -o build_edges_numa
else
	$(CC) -O3 -funroll-loops build_edges.c lplib3.c libmesh6.c -lm -lpthread -o $@
endif

clean:
	rm *.o basic_test basic_test_f77 indirect_writes indirect_writes_f77 tetrahedra_neighbours
