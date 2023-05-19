INCLUDE_DIR=../inc
FLIST=$(basename $(wildcard *.cpp))
OBJS=$(addsuffix .o, $(FLIST))
COMPILER_FLAGS=-g -O0 -Wall -Wextra -pedantic -std=c++17 

# for every objectfile execute recipe
all : $(OBJS)

#Any object-file is created by applying the given command to every cpp file
%.o : %.cpp
	g++ $(COMPILER_FLAGS) -I$(INCLUDE_DIR) -c $< -o $@

#for clean up
clean :
	rm -f *.o

.PHONY: all clean
# DO NOT DELETE

matrix.o: ../inc/matrix.hpp
