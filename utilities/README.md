This directory contains various utilities to convert files in .meshb format to and from some other formats.

### hilbert
This command line renumbers a .mesh(b) file through a Hilbert SFC (Space Filing Curve).
All mesh files entities are renumbered.
This is extremely beneficial to increase cache memory reuse and enhance data paralelism for LPlib multithreading.
