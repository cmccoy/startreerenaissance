CXXFLAGS = -I ../../ -Wall -std=c++11

CPP = $(wildcard src/*.cpp)
OBJ = $(CPP:.cpp=.o)
STATICLIB = libboost-programoptions.a

all: $(STATICLIB)

$(STATICLIB): $(OBJ)
	ar -cvq $@ $^

clean:
	rm -f $(STATICLIB)
	rm -f $(OBJ)

.PHONY: all clean
