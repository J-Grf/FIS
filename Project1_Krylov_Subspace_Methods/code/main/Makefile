INCLUDE_DIR=../inc
EXE=$(addsuffix .exe, $(basename $(wildcard *.cpp)))
COMPILER_FLAGS=-g -O0 -Wall -Wextra -pedantic -std=c++20 

# for every objectfile execute recipe
all : $(EXE)

#Any object-file is created by applying the given command to every cpp file
%.exe : %.cpp
	g++ $(COMPILER_FLAGS) -I$(INCLUDE_DIR) $< ../src/*.o -o $@

#for clean up
clean :
	rm -f $(EXE)

.PHONY: all clean

