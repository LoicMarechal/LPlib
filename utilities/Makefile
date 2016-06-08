

#-----------------------------------------------------------#
#															#
#							HILBERT							#
#															#
#-----------------------------------------------------------#
#															#
#	Description:		makefile LPlib's utilities      	#
#	Author:				Loic MARECHAL						#
#	Creation date:		dec 02 1999							#
#	Last modification:	jun 08 2016							#
#															#
#-----------------------------------------------------------#


# Compiler and flags
CC      = gcc
CFLAGS  = -O3
LDFLAGS = -lm -lpthread

# Exe directories
EXEDIR  = $(HOME)/bin/$(ARCHI)

# Build hilbert
$(EXEDIR)/hilbert: hilbert.c libmesh6.c libmesh6.h lplib3.c lplib3.h
	$(CC) hilbert.c libmesh6.c lplib3.c $(LDFLAGS) $(LDFLAGS) -o $@
