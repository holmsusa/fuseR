
R_HOME = $(R RHOME)
R = $(R_HOME)/bin/R
CC = gcc
#R = R
RM = rm -f

CFLAGS = -g -O2 -fPIC -I$(R_HOME)/include

CFLAGS += -std=c99 -pedantic -Wall -Wextra
# CFLAGS += -march=native
CFLAGS += -fopenmp -fopenmp-simd

CFLAGS += -Wno-unused-function

CFLAGS += -fPIC


LIBS += -lm

## Add flags for profiling with gprof ########
CFLAGS += -pg
LIBS += -pg
##############################################

UNITY_DIR = unity/Unity-master/src
# CFLAGS += -g -pg -fprofile-arcs
# LIBS += -pg -lgcov

all: fuse.so
#testing test
# all: fuse.so test

cleanobj:
	$(RM) fuse.o libfuse.o
#testing.o libfuse.o test.o
# cleanobj:
# 	$(RM) fuse.o libfuse.o test.o

clean: cleanobj
	$(RM) fuse.so
# test testing
# clean: cleanobj
# 	$(RM) fuse.so test

libfuse.o: libfuse.c libfuse.h
	$(CC) $(CFLAGS) -c -o $@ $<

#testing.o: testing.c libfuse.h
#	$(CC) $(CFLAGS) -c -o $@ $< -I$(UNITY_DIR)

#test.o: test.c libfuse.h
#	$(CC) $(CFLAGS) -c -o $@ $<

fuse.so: fuse.c libfuse.o
	PKG_CFLAGS="$(CFLAGS)" PKG_LIBS="$(LIBS)" $(R_HOME)/bin/R CMD SHLIB fuse.c libfuse.o

#testing: testing.o libfuse.o fuse.so
#	$(CC) $(CFLAGS) -o testing testing.o libfuse.o fuse.so -I$(UNITY_DIR) $(UNITY_DIR)/unity.c $(LIBS)

#test: test.o libfuse.o fuse.so
#	$(CC) $(CFLAGS) -o test test.o libfuse.o fuse.so $(LIBS)

.PHONY: all cleanobj clean
.SUFFIXES:
