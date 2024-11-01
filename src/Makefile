
R = R
RM = rm -f

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


# CFLAGS += -g -pg -fprofile-arcs
# LIBS += -pg -lgcov

all: island.so 
# all: island.so test

cleanobj:
	$(RM) island.o libisland.o
# cleanobj:
# 	$(RM) island.o libisland.o test.o

clean: cleanobj
	$(RM) island.so test testing
# clean: cleanobj
# 	$(RM) island.so test

libisland.o: libisland.c libisland.h
	$(CC) $(CFLAGS) -c -o $@ $<

island.so: island.c libisland.o
	PKG_CFLAGS="$(CFLAGS)" PKG_LIBS="$(LIBS)" $(R) CMD SHLIB island.c libisland.o

.PHONY: all cleanobj clean
.SUFFIXES:
