OBJ= newmerge Beam_Ion Beta_gamma TS_check TS_check_raw Beta_ion

all: $(OBJ)

ROOTCFLAGS	=$(shell root-config --cflags)
ROOTLIBS	= $(shell root-config --libs)
ROOTGLIBS	= $(shell root-config --glibs)
TARTLIBS	= -L$(TARTSYS)/lib -lanacore -lanaroot -lanabrips -lanawinds -lanaloop -lanacore
INCLUDES	= -I$(TARTSYS)/include -I$(ROOTSYS)include

CXX	= g++

#BeamGe_AIDA: BeamGe_AIDA.cpp BeamGe_AIDA_mod.h AIDAion.h AIDAbeta.h 
#	$(CXX) $< -o $@ $(INCLUDES) $(TARTLIBS) $(ROOTLIBS) $(ROOTGLIBS) $(ROOTCFLAGS)

newmerge: newmerge.cpp newmerge.h
	$(CXX) $< -o $@ $(INCLUDES) $(TARTLIBS) $(ROOTLIBS) $(ROOTGLIBS) $(ROOTCFLAGS)

Beam_Ion: Beam_ion.cpp
	$(CXX) $< -o $@ $(INCLUDES) $(TARTLIBS) $(ROOTLIBS) $(ROOTGLIBS) $(ROOTCFLAGS)

Beta_ion: Beta_ion.cpp
	$(CXX) $< -o $@ $(INCLUDES) $(TARTLIBS) $(ROOTLIBS) $(ROOTGLIBS) $(ROOTCFLAGS)

Beta_gamma: Beta_gamma.cpp
	$(CXX) $< -o $@ $(INCLUDES) $(TARTLIBS) $(ROOTLIBS) $(ROOTGLIBS) $(ROOTCFLAGS)

TS_check: TS_check.cpp AIDA_ionevt.h AIDA_betaevt.h
	$(CXX) $< -o $@ $(INCLUDES) $(TARTLIBS) $(ROOTLIBS) $(ROOTGLIBS) $(ROOTCFLAGS)

TS_check_raw: TS_check_raw.cpp
	$(CXX) $< -o $@ $(INCLUDES) $(TARTLIBS) $(ROOTLIBS) $(ROOTGLIBS) $(ROOTCFLAGS)

clean:
	rm -f *~ *.o $(OBJ)



