CXXFLAGS = -I ../../ -Wall

CPP = $(wildcard src/*.cpp)
OBJ = $(CPP:.cpp=.o)
STATICLIB = libboost-programoptions.a

all: $(STATICLIB)

$(STATICLIB): $(OBJ)
	ar -cvq $@ $^

.PHONY: all
