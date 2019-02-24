include ../../../Make.global_options
include ../Makefile.options

default : save load

save : save.cpp
	$(MPICXX) -g -o $@ $< $(MPILIBS) $(CXXFLAGS)

./load : load.cpp
	$(MPICXX) -g -o $@ $< $(MPILIBS) $(CXXFLAGS)

clean :
	-rm -rf *~
	-rm -rf *.o
	-rm -rf save load
	-rm *.dx

.PHONY : default clean
