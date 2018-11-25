# Esqueletizacao OpenMP
## Run in serial
./skeletonize ppmimages/elipse.pgm ppmimages/skeleton.pgm 1
## Run in parallel
./skeletonize ppmimages/elipse.pgm ppmimages/skeleton.pgm 0 <nr_threads>
## Compile with papi counters
make papi
## Compile without debugging output
make prod
## Dependencies
gcc 6.5
openmp 4.5
