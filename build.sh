g++ -O3  -std=c++17 -DPARALLEL=0 fista.cpp utilities.cpp -o ./binaries/fastfist
g++ -O3  -std=c++17 greedypp.cpp utilities.cpp -o ./binaries/greedypp
g++ -O3  -std=c++17   mwu.cpp utilities.cpp -o ./binaries/mwu
g++ -O3  -std=c++17   frankwolfe.cpp utilities.cpp -o ./binaries/frankwolfe
g++ -O3  -std=c++17 -DPARALLEL=0 acdm_restart.cpp utilities.cpp -o ./binaries/acdmpr
g++ -O3  -std=c++17 -DPARALLEL=0 rcdm_permutation.cpp utilities.cpp -o ./binaries/rcdmp