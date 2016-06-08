

#-----------------------------------------------------------#
#															#
#							HILBERT							#
#															#
#-----------------------------------------------------------#
#															#
#	Description:		multi-system makefile (gmake only)	#
#	Author:				Loic MARECHAL						#
#	Creation date:		dec 02 1999							#
#	Last modification:	oct 19 2010							#
#															#
#-----------------------------------------------------------#


# default compiler mode is gcc. Other options are: debug, gnu, ibm or hp

MAKEFLAGS = -j 4

# Choose the compiler and its options depending on architecture

ifeq ($(ARCHI), macintel)
	ifeq ($(comp), gnu)
		CC		= /usr/local/bin/gcc
		CFLAGS  = -O3
	endif
	ifeq ($(comp), intel)
		CC		= icc
		CFLAGS  = -fast
		LDFLAGS	= -fast
	endif
	ifeq ($(comp), debug)
        MAKEFLAGS = -j 1
		CC      = gcc
		CFLAGS  = -O -g -Wuninitialized -Wunused
		LDFLAGS = -g
	endif
	ifndef comp
		CC		= /usr/bin/gcc
		CFLAGS  = -O3
		LDFLAGS =
	endif
endif

ifeq ($(ARCHI), macppc)
	ifeq ($(comp), gnu)
		CC		= /usr/local/bin/gcc
		CFLAGS	= -m64 -Ofast
		LDFLAGS	= -m64
	endif
	ifeq ($(comp), ibm)
		CC      = xlc
		CFLAGS  = -O5
		LDFLAGS = -O5
	endif
	ifeq ($(comp), debug)
		CC      = gcc
		CFLAGS  = -O -Wuninitialized -Wunused -g -m64
		LDFLAGS = -g -m64
	endif
	ifndef comp
		CC		= gcc
		CFLAGS  = -fast -m64
		LDFLAGS	= -m64
	endif
endif

ifeq ($(ARCHI), linux64)
	ifeq ($(comp), gnu)
		CC      = gcc
		CFLAGS  = -Ofast -mtune=native -march=native
		LDFLAGS	= -lpthread
	endif
	ifeq ($(comp), intel)
		CC		= icc
		CFLAGS  = -fast -mp1
		LDFLAGS	= -fast -lpthread
	endif
	ifeq ($(comp), debug)
		CC      = gcc
		CFLAGS  = -O -Wuninitialized -Wunused -g
		LDFLAGS = -g -lpthread
	endif
	ifndef comp
		CC      = gcc
		CFLAGS  = -O3 -funroll-loops
		LDFLAGS	= -lpthread
	endif
endif

ifeq ($(ARCHI), linuxppc)
	ifeq ($(comp), gnu)
		CC      = /opt/gcc/bin/gcc
		CFLAGS  = -Ofast
		LDFLAGS	= -lpthread
	endif
	ifeq ($(comp), debug)
		CC      = gcc
		CFLAGS  = -O -Wuninitialized -Wunused -g -m64
		LDFLAGS = -g -m64 -lpthread
	endif
	ifndef comp
		CC      = gcc
		CFLAGS  = -O3 -funroll-loops -m64
		LDFLAGS = -m64 -lpthread
	endif
endif


# Working directories

EXEDIR  = $(HOME)/bin/$(ARCHI)
SRCSDIR = sources
OBJSDIR = objects/$(ARCHI)
ARCHDIR = archives
DIRS    = objects $(EXEDIR) $(OBJSDIR) $(ARCHDIR)
VPATH   = $(SRCSDIR)


# Files to be compiled

SRCS = $(wildcard $(SRCSDIR)/*.c)
HDRS = $(wildcard $(SRCSDIR)/*.h)
OBJS = $(patsubst $(SRCSDIR)%, $(OBJSDIR)%, $(SRCS:.c=.o))
EXE  =  hilbert


# Define compiling implicit rules

$(OBJSDIR)/%.o : $(SRCSDIR)/%.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $(PARFLG) -I$(SRCSDIR) $< -o $@


# Build the main executable

$(EXEDIR)/$(EXE): $(DIRS) $(OBJS)
	$(CC) $(OBJS) $(LDFLAGS) -lm -o $@


# Objects depends on headers

$(OBJS): $(HDRS)


# Build the working directories

$(DIRS):
	@[ -d $@ ] || mkdir $@


# Remove temporary files

clean:
	rm -f $(OBJS)

# Build a dated archive including sources and makefile

tar: $(DIRS)
	tar czf $(ARCHDIR)/$(EXE).`date +"%Y.%m.%d"`.tgz sources Makefile
