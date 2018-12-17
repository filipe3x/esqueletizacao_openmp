# Esqueletizacao OpenMP
## Run in serial
./skeletonize ppmimages/elipse.pgm ppmimages/skeleton.pgm 1
## Run in parallel
./skeletonize ppmimages/elipse.pgm ppmimages/skeleton.pgm 0 <nr_threads>
## Run in parallel with MPI
mpirun -np 4 ./skeletonize ppmimages/horse.pgm ppmimages/skeleton.pgm 4
## Compile with papi counters
make papi
## Compile without debugging output
make prod schedule=<STATIC/DYNAMIC>
## Dependencies
gcc 6.5
openmp 4.5
openmpi 1.8.4
papi 5.4.1
