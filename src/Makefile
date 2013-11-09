CXXFLAGS = -I../lib/eigen_3.2.0 -std=c++0x -Wall -Wextra -g -fopenmp

PROTOBUF_OUT = mutationlist.pb.cc mutationlist.pb.h

ALL = eigen_test $(PROTOBUF_OUT)

all: $(ALL)

gtr.cpp: gtr.hpp

eigen_test: gtr.o mutationlist.pb.o
eigen_test: LDLIBS+= -lprotobuf

%.pb.h %.pb.cc: %.proto
	protoc --cpp_out=. ../python $<

.PHONY: all
