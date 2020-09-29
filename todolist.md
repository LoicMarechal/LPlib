### STANDARD PRIORITY
- develop parallel iterators for FIFO and LIFO stacks
- local scheduling: bind the scheduler to data local to the thread's memory NUMA node
- add a command to kill a pipe while running
- hierarchical block scheduling to enable adaptive block size scheduling
- develop a lattice scheduling based on geometric blocks, not on element indices blocs
- rework the pipeline scheduling so that pipes are queued to a linked list and processed by a small set of threads
- change the pipelines dependency scheme from the current "do not run concurently with pipe X" to "run after completion of pipe X"

### DONE
- handle 64-bit integers
- launchparallel and launchpipeline can take variable numbers of user's arguments
- set and allocate a stack for each thread
