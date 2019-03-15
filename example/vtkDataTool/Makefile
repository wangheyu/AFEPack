# by R.Lie, Nov 01, 2002

include ../../Make.global_options


source = $(wildcard *.cpp)
object = $(patsubst %.cpp, %.o, $(source))
CXXFLAGS = -std=c++11
LDFLAGS += $(LIBS) -ldeal_II
LDFLAGS += -L../../library/lib -lAFEPack

all : main


%.o : %.cpp
	$(CXX) -c -o $@ $< -I./ $(CXXFLAGS)

main : $(object)
	$(CXX) -o $@ $(object) $(LDFLAGS) $(LIBS)
	-@rm -f AFEPack

clean :
	-rm -f $(object)
	-rm -f main
#	-rm -f *.[nes]
#	-@rm -f AFEPack

.PHONY : default clean 

