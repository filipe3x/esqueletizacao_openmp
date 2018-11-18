# Esqueletizacao OpenMP
## Run in serial
./skeletonize ppmimages/elipse.pgm ppmimages/skeleton.pgm 1
## Run in parallel
./skeletonize ppmimages/elipse.pgm ppmimages/skeleton.pgm 0 <nr_threads>
## Dependencies
gcc 6.5
openmp 4.5
