INCLUDE_DIR=../inc
EXE=$(addsuffix .exe, $(basename $(wildcard *.cpp)))
BOOST_LIB = /opt/homebrew/Cellar/boost/1.81.0_1/lib
COMPILER_FLAGS=-g -O0 -Wall -Wextra -pedantic -std=c++17
COMPILER=g++

# for every objectfile execute recipe
all : $(EXE)

#Any object-file is created by applying the given command to every cpp file
%.exe : %.cpp
	$(COMPILER) $(COMPILER_FLAGS) -I$(INCLUDE_DIR) $< ../src/*.o -o $@ -L$(BOOST_LIB) -lboost_timer

#for clean up
clean :
	rm -f *.exe
	rm -f *.txt

.PHONY: all clean

