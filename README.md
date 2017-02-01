#LPlib  version 3.51
A parallelization framework for numerical simulation

#Overview
The purpose of the LPlib is to provide programmers of solvers or automated meshers in the field of scientific computing with an easy, fast and transparent way to parallelize their codes.  
This library is based on posix standard threads, also known as pthreads, thus taking advantage of multi-core chips and shared memory architectures supported by most platforms (Linux, Mac OS X, Unix, Windows).  
It is a simple loop parallelization scheme (hence the name Loop Parallelism Library).  
A serial program can be easily parallelized step by step.  
It requires no knowledge on parallel programing.  
Handles transparently concurrent indirect memory writes and dynamics data structures.

# Build
Simply follow these steps:
-unarchive the ZIP file
-cd LPlib-master
-cmake .
-make
-make install

# Usage
It is made of a single ANSI C file and a header file to be compiled and linked alongside the calling program.  
It may be used in C, C++, Fortran 77 and 90 programs.  
Tested on Linux, Mac OS X, and Windows 7-10.
