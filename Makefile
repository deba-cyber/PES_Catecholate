CXX = g++
CXXFLAGS = -std=c++17 -O3 -malign-double   

#############################

INCDIR = ./include
SRCDIR = ./src
OBJDIR = ./obj
BINDIR = ./bin
OUTPUTDIR = ./output_data

############################

OBJDIRPES := $(OBJDIR)/objpes
OBJDIROPT := $(OBJDIR)/objopt

############################

## Headers  

CONSTHEADER = $(INCDIR)/PES_catecholate/constants.hpp
FILEHEADER = $(INCDIR)/PES_catecholate/fileops.hpp
PESHEADER = $(INCDIR)/PES_catecholate/PES.hpp
OPTHEADER = $(INCDIR)/PES_catecholate/SD.hpp


############################
# Target 

PES: $(OBJDIRPES)/constants.o $(OBJDIRPES)/PES.o $(OBJDIRPES)/main.o
    @echo "Generating PES ..."
    $(CXX) -o $@    $^

$(OBJDIRPES)/constants.o: $(SRCDIR)/constants.cpp $(CONSTHEADER)
    $(CXX) $(CXXFLAGS)  -I$(INCDIR) -c $< -o $@

#$(OBJDIRPES)/PES.o: $(SRCDIR)/PES.cpp $(SRCDIR)/constants.cpp $(CONSTHEADER) $(FILEHEADER) $(PESHEADER)
#   $(CXX) $(CXXFLAGS)  -I$(INCDIR) -c $< -o $@

$(OBJDIRPES)/PES.o: $(SRCDIR)/PES.cpp $(CONSTHEADER) $(FILEHEADER) $(PESHEADER)
    $(CXX) $(CXXFLAGS)  -I$(INCDIR) -c $< -o $@

$(OBJDIRPES)/main.o: $(SRCDIR)/main.cpp $(CONSTHEADER) $(FILEHEADER) $(PESHEADER)
    $(CXX) $(CXXFLAGS)  -I$(INCDIR) -c $< -o $@

#$(OBJDIRPES)/main.o: $(SRCDIR)/main.cpp $(SRCDIR)/PES.cpp $(SRCDIR)/constants.cpp $(CONSTHEADER) $(FILEHEADER) $(PESHEADER)
#   $(CXX) $(CXXFLAGS)  -I$(INCDIR) -c $< -o $@

############################

# Optimisation

OPT: $(OBJDIROPT)/constants.o $(OBJDIROPT)/SD.o $(OBJDIROPT)/PES.o $(OBJDIROPT)/opt.o
    @echo "Finding minimum ..."
    $(CXX) -o $@    $^

$(OBJDIROPT)/constants.o:   $(SRCDIR)/constants.cpp $(CONSTHEADER)
    $(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@

$(OBJDIROPT)/PES.o: $(SRCDIR)/PES.cpp $(CONSTHEADER) $(FILEHEADER) $(PESHEADER)
    $(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@

$(OBJDIROPT)/SD.o: $(SRCDIR)/SD.cpp $(CONSTHEADER) $(OPTHEADER)
    $(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@

$(OBJDIROPT)/opt.o: $(SRCDIR)/opt.cpp $(CONSTHEADER) $(FILEHEADER) $(OPTHEADER) $(PESHEADER)
    $(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@

############################


.PHONY: clean

movepes:
    @echo "Moving PES executable to bin directory ..."
    mv -f PES ./bin

moveopt:
    @echo "Moving OPT executable to bin directory ..."
    mv -f OPT ./bin

clean:
    @echo "Cleaning all exec, objects, binaries ..."
    rm -f $(OBJDIRPES)/*
    rm -f $(OBJDIROPT)/*
    rm -f $(BINDIR)/*
    rm -f $(OUTPUTDIR)/*.bin









