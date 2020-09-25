## LPlib  version 3.62
A parallelization framework for numerical simulation

## Overview
The purpose of the **LPlib** is to provide programmers of solvers or automated meshers in the field of scientific computing with an easy, fast and transparent way to parallelize their codes.  
This library is based on posix standard threads, also known as *pthreads*, thus taking advantage of multi-core chips and shared memory architectures supported by most platforms (*Linux*, *macOS*, *Windows*).  
It is a simple loop parallelization scheme (hence the name Loop Parallelism Library).  
A serial program can be easily parallelized step by step.  
It requires no knowledge on parallel programing.  
Handles transparently concurrent indirect memory writes and dynamics data structures.

## Build

### Prerequisites for *Linux* or *macOS*
- Install [CMake](https://cmake.org/files/v3.7/cmake-3.7.2-win64-x64.msi)
- A valid C99 compiler
- Open a shell window

### Prerequisites for *Windows*
- You first need to install [CMake](https://cmake.org/files/v3.7/cmake-3.7.2-win64-x64.msi). Do not forget to choose "add cmake to the path for all users", from the install panel.
- Then you need a valid C compiler like the free [Visual Studio Community 2019](https://www.visualstudio.com/vs/visual-studio-express/)
- Open the x64 Native Tools Command Prompt for VS (or x86 if you need to build a 32-bit version)

### Build commands for all platforms
- unarchive the ZIP file
- `cd LPlib-master`
- `mkdir build`
- `cd build`
- `cmake ..`
- `cmake --build . --target INSTALL`

### Optional steps
You may download some sample meshes to run the examples:
- you need to install the [libMeshb](https://github.com/LoicMarechal/libMeshb) from GitHub
- manually download files from the *Git LFS* repository: [sample files](sample_meshes/)
- move them into /opt/LPlib/sample_meshes/
- uncompress them with `lzip -d *.meshb.lz`
- you may now enter /opt/LPlib/examples directory and run the various examples

## Usage
It is made of a single *ANSI C* file and a header file to be compiled and linked alongside the calling program.  
It may be used in C, C++, Fortran 77 and 90 programs.  
Tested on *Linux*, *macOS*, and *Windows 7->10*.

Running a parallel loop is pretty easy.  
Let's say that you have a mesh made of vertices and triangles.  
You would like to perform a parallel loop on triangle that would write some data on their three vertices.

Such loop has a _memory race_ since two different threads may process triangles that share a common vertex and write to the same memory location, leading to a wrong result.

Such memory race can be handled by OpenMP with the help of critical sections that you slow down the execution.  
You may also resort to the old "scatter/gather" technic which also slows the execution on top of consuming a lot of extra memory.

Le **LPlib** can handle very efficiently this kind of configuration:

```C++
main()
   // Initialize the library with 4 threads
   LibIdx = InitParallel(4);

   // Set the number of triangles
   TriIdx = NewType(LibIdx, NmbTriangles);

   // Set the number of vertices
   VerIdx = NewType(LibIdx, NmbVertices);

   // Link triangles and their three vertices
   BeginDependency(LibIdx, TriIdx, VerIdx);
   for(i=1;i<=NmbTriangles;i++)
      for(j=0;j<3;j++)
         AddDependency(LibIdx, i, Mesh->Triangles[i][j]);
   EndDependency(LibIdx);

   // Now you can safely launch the parallel loop on triangles
   // telling the library to take care about vertex dependencies
   LaunchParallel(LibIdx, TriIdx, VerIdx, AddSomeValues, Mesh);
}

void AddSomeValues(int begin, int end, int thread, MeshStruct *Mesh)
{
   // Loop over a subset of triangles
   for(i=begin; i<end; i++)
      for(j=0;j<3;j++)
      {
         // Get the vertex index
         VerIdx = Mesh->Triangles[i][j];

         // Modify some vertex' data
         Mesh->vertices[ VerIdx ]->Value += Mesh->Triangle[i]->Value;
         Mesh->vertices[ VerIdx ]->Counter++;
      }
}
```
