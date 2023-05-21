INCLUDE_DIR=../inc
FLIST=matrix #$(basename $(wildcard *.cpp))
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
#gram_schmidt.o: ../inc/gram_schmidt.hpp ../inc/matrix.hpp
#gmres.o: ../inc/gmres.hpp ../inc/matrix.hpp ../inc/gram_schmidt.hpp
