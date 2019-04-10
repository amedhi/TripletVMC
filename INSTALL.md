About
-------
A C++ program for 'Variational Monte Carlo simulation' for quantum lattice 
models.


Version
-------
vmc++ 1.0.0 - Jan 2017.

Build
------------
To build, first make a copy of the file 'make_options.mk' in the 
project root directory and name it as 'options.mk'. Edit this file to set your 
There are two dependencies - Boost and Eigen C++. You will need latest C++ 
compilers as the code contains many C++11 features.  It it tested with 
latest g++ and clang++. 

Usage
-----
To run the application (say, ./a.out), do 

	./a.out [OPTIONS] [FILE]  

The program accepts a few optional arguments OPTIONS. Type 

	./a.out -help

to see the list of available options. The FILE argument, if present, 
is taken as the input filename, else a default filename 'input.parm' is 
looked for in the current directory.  
