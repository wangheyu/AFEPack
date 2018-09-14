include ../Makefile.options

default : save load

save : save.cpp
	$(MPICC) -g -o $@ $< $(MPILIBS)

load : load.cpp
	$(MPICC) -g -o $@ $< $(MPILIBS)

clean :
	-rm -rf *~
	-rm -rf *.o
	-rm -rf save load
	-rm *.dx

.PHONY : default clean
