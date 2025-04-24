
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

UNITY_DIR = unity/Unity-master/src
# CFLAGS += -g -pg -fprofile-arcs
# LIBS += -pg -lgcov

all: island.so 
#testing test
# all: island.so test

cleanobj:
	$(RM) island.o 
#testing.o libisland.o test.o
# cleanobj:
# 	$(RM) island.o libisland.o test.o

clean: cleanobj
	$(RM) island.so
# test testing
# clean: cleanobj
# 	$(RM) island.so test

libisland.o: libisland.c libisland.h
	$(CC) $(CFLAGS) -c -o $@ $<

#testing.o: testing.c libisland.h
#	$(CC) $(CFLAGS) -c -o $@ $< -I$(UNITY_DIR)

#test.o: test.c libisland.h
#	$(CC) $(CFLAGS) -c -o $@ $< 

island.so: island.c libisland.o
	PKG_CFLAGS="$(CFLAGS)" PKG_LIBS="$(LIBS)" $(R) CMD SHLIB island.c libisland.o

#testing: testing.o libisland.o island.so
#	$(CC) $(CFLAGS) -o testing testing.o libisland.o island.so -I$(UNITY_DIR) $(UNITY_DIR)/unity.c $(LIBS)

#test: test.o libisland.o island.so
#	$(CC) $(CFLAGS) -o test test.o libisland.o island.so $(LIBS)

.PHONY: all cleanobj clean
.SUFFIXES:
