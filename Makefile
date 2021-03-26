CXX = g++
CXXFLAGS = -std=c++17 -O3 -malign-double   

#############################

INCDIR = ./include
SRCDIR = ./src
OBJDIR = ./obj
BINDIR = ./bin
OUTPUTDIR = ./output_data

############################

## Headers  

CONSTHEADER = $(INCDIR)/PES_catecholate/constants.hpp 
FILEHEADER = $(INCDIR)/PES_catecholate/fileops.hpp 
PESHEADER = $(INCDIR)/PES_catecholate/PES.hpp


############################

# Target 

PES: $(OBJDIR)/constants.o $(OBJDIR)/PES.o $(OBJDIR)/main.o	
	@echo "Generating PES ..."
	$(CXX) -o $@	$^

$(OBJDIR)/constants.o: $(SRCDIR)/constants.cpp $(CONSTHEADER)
	$(CXX) $(CXXFLAGS)	-I$(INCDIR)	-c $< -o $@

$(OBJDIR)/PES.o: $(SRCDIR)/PES.cpp $(SRCDIR)/constants.cpp $(CONSTHEADER) $(FILEHEADER) $(PESHEADER)
	$(CXX) $(CXXFLAGS)	-I$(INCDIR)	-c $< -o $@

$(OBJDIR)/main.o: $(SRCDIR)/main.cpp $(SRCDIR)/PES.cpp $(SRCDIR)/constants.cpp $(CONSTHEADER) $(FILEHEADER) $(PESHEADER)
	$(CXX) $(CXXFLAGS)	-I$(INCDIR)	-c $< -o $@

############################

.PHONY: clean

move:
	@echo "Moving executables to bin directory ..."
	mv -f PES ./bin

clean:
	@echo "Cleaning all exec, objects, binaries ..."
	rm -f $(OBJDIR)/*
	rm -f $(BINDIR)/*
	rm -f $(OUTPUTDIR)/*.bin





