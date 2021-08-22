# Full Dimensional PES of intramolecular H-transfer in Catecholate

# overview

    PES_Catecholate calculates the potential energy in full dimensional (33)
    normal mode space. The normal modes are dimensionless.  The output is in
    wavenumber (cm^-1)

# Dependencies: C++17 compatibility.

# How to give grid input

    For any given grid, input file is to be generated in ./input_data folder.
    Each line of input file for grid has to have 33 elements corresponding to 33
    dimensionless normal mode displacements.

    In this sample, the file used is ./input_data/catecholate_TS_IRC_grid.dat
    that contains the points along the IRC path

    For each grid, the code needs to be recompiled with the following changes:

    1) In /include/PES_catecholate/constants.hpp ==>

        Set
           inline constexpr int N_pts = Ngrid
        where Ngrid is the no. of grid points for which potential to be
        calculated, i.e. number of lines in the grid input file.

    2) In /src/constants.cpp ==>

        a) Set 
             extern const std::string gridfile = "filename"
           where filename.dat in the ./input_data folder corresponding to grid
           locations. Note: extension .dat is automatically added.

        b) Based on desired output type, following changes are to be made ==>

           If the output energy is desired in binary format, 
           i)  set the following =>
                 extern const std::string outputtype = "bin"
           ii) also set 
                 extern const std::string pot_savefile = "filename"
               where filename.bin will be the output file. Note: extension .bin
               is automatically added.

               After running the executable, output will be saved in the
               ./output_data folder.

           If ascii output is desired:
           i)  set the following
                 extern const std::string outputtype = "dat"
               The potential values will be printed on screen and can be
               redirected to an output file
           ii) also comment out the following line 
                 extern const std::string pot_savefile


# How to run 

    make clean
    make
    make move
    cd bin
    ./PES

   Currently, the code will run for IRC path (80 points) and produce a binary
   output file of energies.

## For any queries
Please contact deba.bhat.90@gmail.com or debabratab@iisc.ac.in
