### Friday the 28th 2025

A new utility to measure memory bandwidth in various mesh-related situations was added.
The command line `cpu\_bandwidth` measure the performance on four different access patterns.

 - Straight access: U( i ) = V( i ) or vector access, the fastest and easiest one to implement.
 - Constant indirect access: U( i ) = T( V( i, 1..4 ) ), memory locations are unpredictable but a constant number of indirect accesses are performed at each step (i.e., accessing the four nodes' data for each tetra)
 - Indirect variable access: U( i ) = T( V(i, 1..d) ), like above but the number of items to be accessed is variable (accessing a vertex's ball of tet in an unstructured mesh)
 - Indirect vectorized variable access: U( i ) = T( V(i, 1..32) ), same as above, but the main loop is split into several sub loops that perform indirect accesses with constant vector sizes with power of two (2,4,8,16,32,64,128 and 256).
 

Command line:
 
```bash
 cpu_bandwidth   NmbThreads   NmbIterations   MeshFile

 In order to fully evaluate all kinds of memory access performances, you should generate
 three different numberings from the same test mesh and feed them to cpu_bandwith:

  hilbert -in MyTestMesh -out RandomMesh -scheme 2
  hilbert -in MyTestMesh -out HilbertMesh
  hilbert -in MyTestMesh -out SlicedMesh -gmlib generic
```
 
How to benchmark your system:
 
Select a tetrahedral test mesh big enough so that cache memory won't spoil the results (bigger than 10 million tets).
Create three different versions as explained above and run the `cpu\_bandwidth` on each mesh.


Interpreting the results:

Random numbering should be fast with direct access and slow with all other kinds of accesses.
Hilbert numbering should be as fast with direct access and almost as fast with indirect accesses. It means that with an efficient renumbering, all indirect accesses created by unstructured meshes are almost as fast as direct access used by structures meshes.
Variable indirect accesses are a little slower than constant indirect accesses and even vectorized variable accesses are not any faster on CPUs. Such a method is only efficient with GPUs.
Parallelism: direct memory access scales up to the number of memory buses and stalls above while indirect accesses scale quite well.
