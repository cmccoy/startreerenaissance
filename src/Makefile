CXXFLAGS = -I../lib/eigen_3.2.0 -std=c++0x -Wall -Wextra -g

PROTOBUF_OUT = mutationlist.pb.cc mutationlist.pb.h

EXE = fit_gtr build_mutation_matrices

ALL = $(EXE) $(PROTOBUF_OUT)

all: $(ALL)

clean:
	rm -f $(EXE) *.o

gtr.o: gtr.hpp

fit_gtr: gtr.o mutationlist.pb.o
fit_gtr: LDLIBS += -lprotobuf -lz
fit_gtr: CXXFLAGS += -fopenmp

# Protobuf stuff
%.pb.h %.pb.cc: %.proto
	protoc --cpp_out=. ../python $<

# BAM
LIBBAM = ../lib/samtools/libbam.a
$(LIBBAM):
	make -C ../lib/samtools libbam.a

# Mutation matrices
build_mutation_matrices: mutationlist.pb.o $(LIBBAM)
build_mutation_matrices: LDLIBS += -lz -lpthread -lprotobuf
build_mutation_matrices: CXXFLAGS += -I../lib/samtools

.PHONY: all clean



# Gross auto-dependencies
#%.o : %.cc
	#@g++ -MD -c -o $@ $<
	#@cp $*.d $*.P; \
			 #sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
					 #-e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.P; \
			 #rm -f $*.d

#-include *.P
