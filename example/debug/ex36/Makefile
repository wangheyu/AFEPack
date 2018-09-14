include ../Makefile.options

default : trilinos petsc

trilinos : trilinos.cpp
	$(MPICC) -o $@ $< $(CFLAGS) $(MPILIBS) $(TRILINOS_INC) $(TRILINOS_LIBS)

petsc : petsc.cpp
	$(PETSC_MPICC) -o $@ $< $(CFLAGS) $(MPILIBS) $(PETSC_INC) $(PETSC_LIBS)

clean :
	-rm -rf *~
	-rm -rf *.o
	-rm -rf trilinos petsc
	-rm *.dx

.PHONY : default clean
