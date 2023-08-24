CXX     = g++

CXXFLAGS= $(shell root-config --cflags)
LIBS    = $(shell root-config --libs) 


SOURCES = ./Main_calib.cc
HEADERS = ./Main_calib.h
OBJECTS = $(SOURCES:.cc=.o)

EXECUTABLE = PFCalib

all: $(SOURCES) $(EXECUTABLE)

%.o: %.cc
	@echo Compiling $<...
	$(CXX) $(CXXFLAGS) -c -o $@ $< 


$(EXECUTABLE): $(OBJECTS)
	@echo "Linking $(PROGRAM) ..."
	@echo "@$(CXX) $(LIBS) $(OBJECTS) -o $@"
	@$(CXX) -o $@ $^ $(LIBS)
	@echo "done"


# Specifying the object files as intermediates deletes them automatically after the build process.
.INTERMEDIATE:  $(OBJECTS)

# The default target, which gives instructions, can be called regardless of whether or not files need to be updated.
.PHONY : clean
clean:
	rm -f $(OBJECTS) $(EXECUTABLE)

###
Main_calib.o: ./Main_calib.h
# RecHit.o: ./interface/RecHit.h
# LateralRings.o: ./interface/LateralRings.h
# RecHit_AHCAL.o: ./interface/RecHit_AHCAL.h
# LateralSquares.o: ./interface/LateralSquares.h
# HGCNtupleVariables.o: ./interface/HGCNtupleVariables.h
# AnalyzeHGCOctTB.o: ./interface/RecHit.h ./interface/LateralRings.h ./interface/HGCNtupleVariables.h ./interface/AnalyzeHGCOctTB.h
