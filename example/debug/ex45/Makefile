include ../../../Make.global_options
include ../Makefile.options

source = $(wildcard *.cpp)
object = $(patsubst %.cpp, %.o, $(source))

default : main

%.o : %.cpp
	$(MPICXX) -c -o $@ $< $(CXXFLAGS)

main : $(object)
	$(MPICXX) -o $@ $(object) $(MPILIBS)

clean :
	-rm -rf run.sh
	-rm -rf *~
	-rm -rf $(object)
	-rm -rf main
	-rm *.dx

.PHONY : default clean
