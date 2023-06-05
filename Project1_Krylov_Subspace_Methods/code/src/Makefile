INCLUDE_DIR = ../inc
BOOST_PATH = /opt/homebrew/Cellar/boost/1.81.0_1/include/
FLIST = $(basename $(wildcard *.cpp))
OBJS = $(addsuffix .o, $(FLIST))
COMPILER_FLAGS = -O0 -Wall -Wextra -pedantic -std=c++17 #-DDISABLEIO
COMPILER = g++

# for every objectfile execute recipe
all : $(OBJS)

#Any object-file is created by applying the given command to every cpp file
%.o : %.cpp
	$(COMPILER) $(COMPILER_FLAGS) -I$(INCLUDE_DIR) -I$(BOOST_PATH) -c $< -o $@

#for clean up
clean :
	rm -f *.o

.PHONY: all clean
