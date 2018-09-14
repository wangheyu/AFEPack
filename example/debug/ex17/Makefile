include ../Makefile.options

source = $(wildcard *.cpp)
object = $(patsubst %.cpp, %.o, $(source))

default : main

%.o : %.cpp
	$(MPICC) -c -o $@ $< $(CFLAGS)

main : $(object)
	$(MPICC) -o $@ $(object) $(MPILIBS)

clean :
	-rm -rf run.sh
	-rm -rf *~
	-rm -rf $(object)
	-rm -rf main
	-rm *.dx

.PHONY : default clean
