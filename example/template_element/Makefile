# by R.Lie, Nov 01, 2002

# enable debug under default case
ifeq (x$(debug), x)
  debug=yes
endif

source = $(wildcard *.cpp)
ifeq ($(debug), yes)
  object = $(patsubst %.cpp, %.g.o, $(source))
  CXXFLAGS = -g -DDEBUG
  lib = -lAFEPack.g -llac.g -lbase.g -ldl -lpthread $(EXTRA_LIBDIR)
else
  object = $(patsubst %.cpp, %.o, $(source))
  CXXFLAGS = -O6
  lib = -lAFEPack -llac -lbase -ldl -lpthread $(EXTRA_LIBDIR)
endif

all : main

%.g.o : %.cpp
	$(CXX) -c -o $@ $< $(CXXFLAGS) $(EXTRA_INCDIR)

%.o : %.cpp
	$(CXX) -c -o $@ $< $(CXXFLAGS) $(EXTRA_INCDIR)

main : $(object)
	$(CXX) -o $@ $(object) $(LDFLAGS) $(lib)

clean :
	-rm -rf $(object)
	-rm -rf main
	-rm -f *.[nes]
	-rm -f *.dx

.PHONY : default clean
