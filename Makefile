##########################################################################
# Makefile for sparse coordinate
##########################################################################

HDR = simple_sparse_vec_hash.h WeightVector.h coordinate_descent.h
SRC = simple_sparse_vec_hash.cpp  WeightVector.cpp coordinate_descent.cpp cmd_line.cpp

CC = g++

CC      = g++
CFLAGS  = -Wall -O3 
#CFLAGS  = -g 
LFLAGS  = -lm

OBJS = $(SRC:.cpp=.o)

all: coordinate test_objective

test_objective: $(OBJS) test_objective.o
	$(CC) $(OBJS) test_objective.o $(LFLAGS) -o test_objective

coordinate: $(OBJS) main.o
	$(CC) $(OBJS) main.o $(LFLAGS) -o coordinate

tar: $(SRC) $(HDR) Makefile 
	tar zcvf coordinate.tgz *.cpp *.h Makefile license.txt README data.tgz

%.o: %.cpp $(HDR)
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o *.od *.oc *~ \#*\# depend coordinate coordinate.exe test_objective.exe
