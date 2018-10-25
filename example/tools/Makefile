# by R.Lie, Nov 01, 2002

include ../../Make.global_options

install-dir = $(prefix)
lib = -L../../library/lib -lAFEPack $(LDFLAGS) $(LIBS)

with_script = \
  2d_vecfun_dxdata2tecplotdata.bin \
  dbmesh2mesh.bin \
  dxdata2easymeshdata.bin \
  dxdata2tecplotdata.bin \
  gmsh2mesh.bin 

indenped = \
  easymesh2mesh.exe \
  easymesh_refine.exe \
  easymesh2opendx.exe \
  mesh2easymesh.exe \
  mesh2opendx.exe \
  mesh_refine.exe \
  hsfc_renum.exe \
  tecplotdata2dxdata.exe

executable = $(with_script) $(indenped)

all : $(executable)

%.bin : %.cpp
	$(CXX) -o $@ $< $(CXXFLAGS) $(LDFLAGS) $(lib)

%.exe : %.cpp
	$(CXX) -o $@ $< $(CXXFLAGS) $(LDFLAGS) $(lib)
	@mv $@ `echo $@ | sed -e "s/.exe//g"`

install:
	@list=`ls`; for file in $$list; do \
	  mkdir -p $(install-dir)/bin ; \
	  if test -x $$file; then \
	    cp $$file $(install-dir)/bin ; \
	  fi; \
	done

clean :
	-rm -f $(patsubst %.exe, %, $(executable))
	-rm -f *~

.PHONY : default clean
