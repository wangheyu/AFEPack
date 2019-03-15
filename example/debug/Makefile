include ../../Make.global_options
include Makefile.options

# ex17 ex19, ex20, ex31 ex32 ex33: 守恒律有限体积并行例子
# ex24 ex25: Trilinos 例子
# ex26 ex27: PETSc 例子
# ex28 ex29 ex30: 背景网格分区
# ex38 ex39: HYPRE 例子
# ex42 ex43 ex44: 可以调整根单元的负载平衡的例子
# ex45: ex20改为可以调整根单元的负载平衡的情况

EXE=ex03 ex04 

MPIEXE=ex05 ex06 ex07 ex08 ex09 ex10 \
       ex11 ex12 ex13 ex14 ex15 ex16 ex18 ex19 \
       ex28 ex29-2d ex29-3d ex30 ex41 \
       ex42 ex43 ex44

TRILIONS=ex24 ex25

PETSC=ex26 ex27

HYPRE=ex38 ex39

SUBDIRS=ex17 ex19 ex20 ex31 ex32 ex33 ex34 ex35 ex36 ex37 ex40 ex45

all : $(EXE) $(MPIEXE) $(TRILIONS) $(PTESC)
	@list='$(SUBDIRS)'; for i in $$list ; do \
	  ( cd $$i && make ) \
	done

ex03 : ex03.cpp
	$(MPICXX) -o $@ $< $(CXXFLAGS) $(LIBS)

ex04 : ex04.cpp
	$(MPICXX) -o $@ $< $(CXXFLAGS) $(LIBS)

ex05 : ex05.cpp
	$(MPICXX) -o $@ $< $(CXXFLAGS) $(MPILIBS)

ex06 : ex06.cpp
	$(MPICXX) -o $@ $< $(CXXFLAGS) $(MPILIBS)

ex07 : ex07.cpp
	$(MPICXX) -o $@ $< $(CXXFLAGS) $(MPILIBS)

ex08 : ex08.cpp
	$(MPICXX) -o $@ $< $(CXXFLAGS) $(MPILIBS)

ex09 : ex09.cpp
	$(MPICXX) -o $@ $< $(CXXFLAGS) $(MPILIBS)

ex10 : ex10.cpp
	$(MPICXX) -o $@ $< $(CXXFLAGS) $(MPILIBS)

ex11 : ex11.cpp
	$(MPICXX) -o $@ $< $(CXXFLAGS) $(MPILIBS)

ex12 : ex12.cpp
	$(MPICXX) -o $@ $< $(CXXFLAGS) $(MPILIBS)

ex13 : ex13.cpp
	$(MPICXX) -o $@ $< $(CXXFLAGS) $(MPILIBS)

ex14 : ex14.cpp
	$(MPICXX) -o $@ $< $(CXXFLAGS) $(MPILIBS)

ex15 : ex15.cpp
	$(MPICXX) -o $@ $< $(CXXFLAGS) $(MPILIBS)

ex16 : ex16.cpp
	$(MPICXX) -o $@ $< $(CXXFLAGS) $(MPILIBS)

ex17 : 
	cd $@ && make

ex18 : ex18.cpp
	$(MPICXX) -o $@ $< $(CXXFLAGS) $(MPILIBS)

ex19 :
	cd $@ && make

ex20 :
	cd $@ && make

ex24 : ex24.cpp
	$(MPICXX) -o $@ $< $(CXXFLAGS) $(MPILIBS) $(TRILINOS_INC) $(TRILINOS_LIBS)

ex25 : ex25.cpp
	$(MPICXX) -o $@ $< $(CXXFLAGS) $(MPILIBS) $(TRILINOS_INC) $(TRILINOS_LIBS)

ex26 : ex26.cpp
	$(PETSC_MPICXX) -o $@ $< $(CXXFLAGS) $(MPILIBS) $(PETSC_INC) $(PETSC_LIBS)

ex27 : ex27.cpp
	$(PETSC_MPICXX) -o $@ $< $(CXXFLAGS) $(MPILIBS) $(PETSC_INC) $(PETSC_LIBS)

ex28 : ex28.cpp
	$(MPICXX) -o $@ $< $(CXXFLAGS) $(MPILIBS)

ex29-2d : ex29.cpp
	$(MPICXX) -o $@ $< $(CXXFLAGS) $(MPILIBS)

ex29-3d : ex29.cpp
	$(MPICXX) -DDIM3 -o $@ $< $(CXXFLAGS) $(MPILIBS)

ex30 : ex30.cpp
	$(MPICXX) -o $@ $< $(CXXFLAGS) $(MPILIBS)

ex38 : ex38.cpp
	$(HYPRE_MPICXX) -o $@ $< -fPIC $(CXXFLAGS) $(HYPRE_LIBS) $(HYPRE_INC) $(MPILIBS) 

ex39 : ex39.cpp
	$(HYPRE_MPICXX) -o $@ $< -fPIC $(CXXFLAGS) $(HYPRE_LIBS) $(HYPRE_INC) $(MPILIBS) 

ex41 : ex41.cpp
	$(MPICXX) -o $@ $< $(CXXFLAGS) $(MPILIBS)

ex42 : ex42.cpp
	$(MPICXX) -o $@ $< $(CXXFLAGS) $(MPILIBS)

ex43 : ex43.cpp
	$(MPICXX) -o $@ $< $(CXXFLAGS) $(MPILIBS)

ex44 : ex44.cpp
	$(MPICXX) -o $@ $< $(CXXFLAGS) $(MPILIBS)

clean :
	-rm $(EXE)
	-rm $(MPIEXE)
	-rm $(TRILIONS)
	-rm $(PETSC)
	-rm $(HYPRE)
	-rm *~
	-rm *.dx
	@list='$(SUBDIRS)'; for i in $$list ; do \
	  ( cd $$i && make clean ) \
	done

.PHONY : clean $(SUBDIRS)

