INCLUDE_DIR = ../inc
FLIST = $(basename $(wildcard *.cpp))
OBJS = $(addsuffix .o, $(FLIST))
COMPILER_FLAGS =-march=native -O3 -Wall -Wextra -pedantic -std=c++20 -DDISABLEIO
COMPILER = g++

# for every objectfile execute recipe
all : $(OBJS)

#Any object-file is created by applying the given command to every cpp file
%.o : %.cpp
	$(COMPILER) $(COMPILER_FLAGS) -I$(INCLUDE_DIR) -c $< -o $@

#for clean up
clean :
	rm -f *.o

.PHONY: all clean
