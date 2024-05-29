
Implementations for Ta Duy Nguyen & Alina Ene:

**Multiplicative Weights Update, Area Convexity and Random Coordinate Descent for Densest Subgraph Problems**

Codebase is credited to Harb, Quanrud, Chekuri (NeurIPS 2022): **Faster and Scalable Algorithms for Densest Subgraph and Decomposition**

## General structure & Compilation Instructions

To compile, use build.sh

## Input format
The input format is expected to be a .txt file, of the format 
n m
u1 v1
u2 v2
...
um vm

where n is the number of vertices in G, m is the number of edges, and 0<= ui,vi <= n-1 for all i. Ensure there are no self loops (the algorithms were not tested with self loops, they might work, but care needed). 

In addition, if the optimum load vector is known (i.e b^\ast is known), then the following format is allowed:
n m
u1 v1
u2 v2
...
um vm
b0
b1
...
bn-1

There bi is the i-th vertex optimum density. The algorithm would output the error with respect to bi given at every iteration. 


## Example run 

After compiling, you can test the algorithm on the "Close cliques" dataset using the following:

$ ./binaries/fastfist 200 < ./datasets_cleaned/close-cliques.txt


The 200 here signifies the number of iterations to run the algorithm. All binaries support the same syntax (for example change fastfist above with mwu and it still works using the MWU algorithm). 
