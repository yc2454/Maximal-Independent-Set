# Maximal-Independent-Set

This repo contains C++ implementations of algorithms for finding Maximal Independent Sets, whose definition can be found here: https://en.wikipedia.org/wiki/Maximal_independent_set.

`MIS_seq.cpp` contains the implementation of a simple serial algorithm for finding all MIS in a graph.

`MIS_para.cpp` contains the implementation of Luby's Algorithm, a parallel algorithm with a span of O(log n). The implementation makes use of **Open MPI**, version 4.0.5, an open source parallel computing library for C++. Details about the library can be found here: https://www.open-mpi.org/.
