# Parallel-Graph-Bisection
A parallel implementation of the label propagation-based graph bisection algorithm using MPI! Third year at ENS project done in 2020

# How to launch
  - Make sure to have MPI installed. 
  - In the "/code" repository execute the "make" commande to compile the project.  
  - Now execute the program as follow:   mpirun -np [num-procs] ./graph-bisect [graph-file-name] [bisect-algo-name] [max-iters] [epsilon] [random-seed]    
        - [num-procs]: The number of MPI ranks to be used. For the sequential algorithm it has to be 1.  
        - [graph-file-name]: Name of the graph file.  
        - [bisect-algo-name]: Name of the graph bisection algorithm. There are three possibilities:  
          * bisect-seq: Sequential bisection algorithm.  
          * bisect-a2a: Parallel bisection algorithm using all-to-all communication routines.  
          * bisect-p2p: Parallel bisection algorithm using point-to-point communication routines.  
        - [max-iters]: The number of iterations for the label propagation to run.  
        - [epsilon]: Maximum allowed imbalance ratio for the label propagation.  
        - [seed]: Random seed used to initialize the part array. This parameter is optional, and you can use it to test the correctness of your algorithm. If not specified, it is set to time(0).  
        
   
