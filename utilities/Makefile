

#-----------------------------------------------------------------------------#
#                                                                             #
#                                LPlib utilities                              #
#                                                                             #
#-----------------------------------------------------------------------------#
#                                                                             #
#  Description:         makefile LPlib's utilities                            #
#  Author:              Loic MARECHAL                                         #
#  Creation date:       mar 02 2017                                           #
#  Last modification:   jul 08 2020                                           #
#                                                                             #
#-----------------------------------------------------------------------------#


# Compiler and flags
CC      = clang
CFLAGS  = -Ofast -flto -I ~/include -I . -DWITH_GMLIB
LDFLAGS = -lm

# Exe directories
EXEDIR  = $(HOME)/bin/$(ARCHI)

# Build hilbert
$(EXEDIR)/hilbert: hilbert.c libmeshb7.c libmeshb7.h lplib3.c lplib3.h
	$(CC) hilbert.c libmeshb7.c lplib3.c $(CFLAGS) $(LDFLAGS) -o $@

clean:
	rm $(EXEDIR)/hilbert
