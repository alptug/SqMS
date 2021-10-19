#
# 'make depend' uses makedepend to automatically generate dependencies 
#               (dependencies are added to end of Makefile)
# 'make'        build executable file 'mycc'
# 'make clean'  removes all .o and executable files
#

# define the C compiler to use
CC=gcc

# define any compile-time flags
CFLAGS=-Wall -O0  -g -pg

# define any directories containing header files other than /usr/include
#
INCLUDES=#-I./include

# define library paths in addition to /usr/lib
#   if I wanted to include libraries not in /usr/lib I'd specify
#   their path using -Lpath, something like:
LFLAGS=#-L./lib

# define any libraries to link into executable:
#   if I want to link in libraries (libx.so or libx.a) I use the -llibname 
#   option, something like (this will link in libmylib.so and libm.so:
LIBS=-lm

# define the source code directory
SRC=src
# define the object directory
OBJ=obj

# find the C source files
SRCS=$(wildcard $(SRC)/*.c)
# create objects
OBJS=$(patsubst $(SRC)/%.c, $(OBJ)/%.o, $(SRCS))


# define the executable file 
BINDIR=bin
BIN = $(BINDIR)/main

DATADIR=data
#
# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#

.PHONY: all clean

all: $(BIN)
	@echo Project has been compiled

release: CFLAGS=-Wall -O3 -DNDEBUG -march=native
release: clean
release: $(BIN)

$(BIN): $(OBJS) 
	$(CC) $(CFLAGS) $(INCLUDES) $(OBJS) -o $@  $(LFLAGS) $(LIBS)

# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file) 
# (see the gnu make manual section about automatic variables)
$(OBJ)/%.o: $(SRC)/%.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	$(RM) -r $(BINDIR)/* $(OBJ)/* $(DATADIR)/*

#depend: $(SRCS)
#		makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it