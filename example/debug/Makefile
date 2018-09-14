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
	$(CC) -o $@ $< $(CFLAGS) $(LIBS)

ex04 : ex04.cpp
	$(CC) -o $@ $< $(CFLAGS) $(LIBS)

ex05 : ex05.cpp
	$(MPICC) -o $@ $< $(CFLAGS) $(MPILIBS)

ex06 : ex06.cpp
	$(MPICC) -o $@ $< $(CFLAGS) $(MPILIBS)

ex07 : ex07.cpp
	$(MPICC) -o $@ $< $(CFLAGS) $(MPILIBS)

ex08 : ex08.cpp
	$(MPICC) -o $@ $< $(CFLAGS) $(MPILIBS)

ex09 : ex09.cpp
	$(MPICC) -o $@ $< $(CFLAGS) $(MPILIBS)

ex10 : ex10.cpp
	$(MPICC) -o $@ $< $(CFLAGS) $(MPILIBS)

ex11 : ex11.cpp
	$(MPICC) -o $@ $< $(CFLAGS) $(MPILIBS)

ex12 : ex12.cpp
	$(MPICC) -o $@ $< $(CFLAGS) $(MPILIBS)

ex13 : ex13.cpp
	$(MPICC) -o $@ $< $(CFLAGS) $(MPILIBS)

ex14 : ex14.cpp
	$(MPICC) -o $@ $< $(CFLAGS) $(MPILIBS)

ex15 : ex15.cpp
	$(MPICC) -o $@ $< $(CFLAGS) $(MPILIBS)

ex16 : ex16.cpp
	$(MPICC) -o $@ $< $(CFLAGS) $(MPILIBS)

ex17 : 
	cd $@ && make

ex18 : ex18.cpp
	$(MPICC) -o $@ $< $(CFLAGS) $(MPILIBS)

ex19 :
	cd $@ && make

ex20 :
	cd $@ && make

ex24 : ex24.cpp
	$(MPICC) -o $@ $< $(CFLAGS) $(MPILIBS) $(TRILINOS_INC) $(TRILINOS_LIBS)

ex25 : ex25.cpp
	$(MPICC) -o $@ $< $(CFLAGS) $(MPILIBS) $(TRILINOS_INC) $(TRILINOS_LIBS)

ex26 : ex26.cpp
	$(PETSC_MPICC) -o $@ $< $(CFLAGS) $(MPILIBS) $(PETSC_INC) $(PETSC_LIBS)

ex27 : ex27.cpp
	$(PETSC_MPICC) -o $@ $< $(CFLAGS) $(MPILIBS) $(PETSC_INC) $(PETSC_LIBS)

ex28 : ex28.cpp
	$(MPICC) -o $@ $< $(CFLAGS) $(MPILIBS)

ex29-2d : ex29.cpp
	$(MPICC) -o $@ $< $(CFLAGS) $(MPILIBS)

ex29-3d : ex29.cpp
	$(MPICC) -DDIM3 -o $@ $< $(CFLAGS) $(MPILIBS)

ex30 : ex30.cpp
	$(MPICC) -o $@ $< $(CFLAGS) $(MPILIBS)

ex38 : ex38.cpp
	$(HYPRE_MPICC) -o $@ $< -fPIC $(CFLAGS) $(HYPRE_LIBS) $(MPILIBS) 

ex39 : ex39.cpp
	$(HYPRE_MPICC) -o $@ $< -fPIC $(CFLAGS) $(HYPRE_LIBS) $(MPILIBS) 

ex41 : ex41.cpp
	$(MPICC) -o $@ $< $(CFLAGS) $(MPILIBS)

ex42 : ex42.cpp
	$(MPICC) -o $@ $< $(CFLAGS) $(MPILIBS)

ex43 : ex43.cpp
	$(MPICC) -o $@ $< $(CFLAGS) $(MPILIBS)

ex44 : ex44.cpp
	$(MPICC) -o $@ $< $(CFLAGS) $(MPILIBS)

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

