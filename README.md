#LPlib  version 3.51
A parallelization framework for numerical simulation

#Overview
The purpose of the **LPlib** is to provide programmers of solvers or automated meshers in the field of scientific computing with an easy, fast and transparent way to parallelize their codes.  
This library is based on posix standard threads, also known as *pthreads*, thus taking advantage of multi-core chips and shared memory architectures supported by most platforms (*Linux*, *macOS*, *Windows*).  
It is a simple loop parallelization scheme (hence the name Loop Parallelism Library).  
A serial program can be easily parallelized step by step.  
It requires no knowledge on parallel programing.  
Handles transparently concurrent indirect memory writes and dynamics data structures.

# Build
Simply follow these steps:
- unarchive the ZIP file
- `cd LPlib-master`
- `cmake .`
- `make`
- `make install`

Optionally, you may download some sample meshes to run the examples:
- you need to install the [libMeshb](https://github.com/LoicMarechal/libMeshb) from GitHub
- manually download files from the *Git LFS* repository: [sample files](sample_meshes/)
- move them into /opt/LPlib/sample_meshes/
- uncompress them with `lzip -d *.meshb.lz`
- you may now enter /opt/LPlib/examples directory and run the various examples

# Usage
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
   for(i=1;i<=NmbTriangles;i++)
      for(j=0;j<3;j++)
         AddDependency(LibIdx, i, Mesh->Triangles[i][j])
   EndDependency(LibIdx);

   // Now you can safely launch the parallel loop on triangles
   // telling the library to take care about vertex dependencies
   LaunchParallel(LibIdx, TriIdx, VerIdx, AddSomeValue, Mesh);
}

void AddSomeValue(int begin, int end, int thread, MeshStruct *Mesh)
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
