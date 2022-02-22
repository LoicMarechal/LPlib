### STANDARD PRIORITY
- develop parallel iterators for FIFO and LIFO stacks
- local scheduling: bind the scheduler to data local to the thread's memory NUMA node
- add a command to kill a pipe while running
- hierarchical block scheduling to enable adaptive block size scheduling
- develop a lattice scheduling based on geometric blocks, not on element indices blocs
- link dependency block at creation and do not unlink them while running the parallel loop
- interleaved procedures: allow multiple procedures to be launched in parallel and processed in a pipelined way

### DONE
- handle 64-bit integers
- launchparallel and launchpipeline can take a variable number of user's arguments
- set and allocate a stack for each thread
- set and allocate a stack for each pipe
- change the pipeline dependency scheme from the current "do not run concurrently with pipe X" to "run after completion of pipe X"
- interleave independent data blocks to mitigate threads unbalance
- rework the pipeline scheduling so that pipes are queued to a linked list and processed by a small set of threads
- add an option to choose between a fixed number of blocks and a fixed block size in interleaved mode
- make interleaved blocks optional as some parallel loops cannot work with it (i.e. TetrahedraNeighbours example)
- add a procedure to set new parameters before calling LaunchParallel()
- develop an optional static scheduling option that makes the process deterministic
