# Add source files
EXECUTABLE := 1dwavesolver

CFILES := 1dwavesolver.cpp
CDEPS  :=
HFILES := 
CFLAGS := -lgsl -lgslcblas -lfftw3 -std=c++17 -lm 

# Rules and targets
$(EXECUTABLE) : $(CFILES) $(CDEPS) $(HFILES)
	g++ -o $@ $(CFILES) $(CFLAGS)