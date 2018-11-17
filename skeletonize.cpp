#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <papi.h>

#define BLACKPIXEL 0
#define WHITEPIXEL 1

/*

Call:
void skeletonize(img->buf, img->width, img->height) 

The strategy is to aglomerate data and therefore, computations. Ideally we should have each data chunck fitting L2 cache. Each chunck should be given to a thread/core.
We can check if memory calls are flooding to L3 with papi counters.

*/

void skeletonize_serial(int *I, int W, int H) {
	int *neighbors = (int*) malloc(9*sizeof(int)); // each pixel will have 8 neighbors
	int *chan1to0 = (int*) malloc(W*H*sizeof(int)); // which pixels we are going to change
	int *cont = (int*) malloc(2*sizeof(int)); // for checking if we already finished

	int X_index[8] = {-1,-1,0,1,1,1,0,-1}; // neighbors relative coordinates
	int Y_index[8] = {0,1,1,1,0,-1,-1,-1};

	int i, j, k; // indexes
	int total; // total of neighbors
	int ans; // total transitions from 0 to 1

	long long PAPI_start, PAPI_stop; 
	double stop, start = omp_get_wtime();
	PAPI_start = PAPI_get_real_usec();

	// each thread will have its own value of x, total, i, j
	// #pragma omp parallel for private(x, total, i, j) schedule(static)
	while(cont[0] > 0 || cont[1] > 0) {
		cont[0] = 0;
		cont[1] = 0;

		for(i=1; i < H-1; i++) {
			for(j=1; j < W-1; j++) {
				total = 0;
				ans = 0;
				chan1to0[i*H+j] = 0;

				for(k=0; k < 8; k++) { // get all neighbors
					neighbors[k] = I[(i+X_index[k])*W + (j+Y_index[k])];
					total += neighbors[k];
				}
				neighbors[8] = total;

				for(k=0; k < 7; k++) { // count nr transitions from 0 to 1
					if(neighbors[k] == 0 && neighbors[k+1] == 1)
						ans += 1;
				}

				if(neighbors[7] == 0 && neighbors[0] == 1) ans += 1;

				if(i % 2 != 0 && I[i*W+j] == 1 && neighbors[8] >= 2 && neighbors[8] <= 6 && ans == 1 
					&& neighbors[0] * neighbors[2] * neighbors[4] == 0 
					&& neighbors[2] * neighbors[4] * neighbors[6] == 0) {
					
					chan1to0[i*H+j] = 1; // we mark pixel for deletion
					cont[0] = 1;

				}

				if(i % 2 == 0 && I[i*W+j] == 1 && neighbors[8] >= 2 && neighbors[8] <= 6 && ans == 1 
					&& neighbors[0] * neighbors[2] * neighbors[6] == 0 
					&& neighbors[0] * neighbors[4] * neighbors[6] == 0) {
					
					chan1to0[i*H+j] = 1; // we mark pixel for deletion
					cont[1] = 1;

				}
			}
		}

		// we delete pixels
		for(i=1; i < H-1; i++) {
			for(j=1; j < W-1; j++) {
				if(chan1to0[i*W+j] == 1) I[i*W+j] = 0;
			}
		}
	}
	stop = omp_get_wtime();
	PAPI_stop = PAPI_get_real_usec();
	printf("papi Time in microseconds: %lld\n", PAPI_stop - PAPI_start);
	printf("omp Time in microseconds: %f\n",(stop - start)*1000000 );
}

/*
// i -> line, j -> collum
int* posicao(int i, int j, int** img) {
	return &img[i][j];	
}

int *getNeighbors(int *pixel, int p, int W) {
	int* neighbors = malloc(sizeof(int)*8); // each pixel has 8 neighbors, for a mask of size 3
	neighbors[0] = pixel[p - W*p]; // P2
	neighbors[1] = pixel[p + 1]; // P3
	neighbors[2] = pixel[p + W*p +1]; // P4
	neighbors[3] = pixel[p + W*p]; // P5
	neighbors[4] = *posicao(i + 1,p); // P6
	neighbors[5] = *posicao(i + 1,p - 1); // P7
	neighbors[6] = *posicao(i, p - 1); // P8
	neighbors[7] = *posicao(i - 1,p - 1); // P9
	return neighbors;	
}

int neighborsNumber(int* neighbors) {
	int ones = 0;
	for(int i = 0; i <= 8; i++)
		neighbors[i] ? ones++ : 0; 	
	return ones;
}

int sequenceCounter(int* neighbors) {
	int count = 0;
	int i = 0;
	for( ; i < 8; i++)
		if(neighbors[i] < neighbors[i+1]) count++;

	if(neighbors[i+1] < neighbors[0]) count++;

	return count;
}

int firstPassDiagonal(int* neighbors) {
	return neighbors[4] * neighbors[6] * (neighbors[0] + neighbors[2])
}

int secondPassDiagonal(int* neighbors) {
	return neighbors[0] * neighbors[2] * (neighbors[4] + neighbors[6])
}

int addPadding(int* img, int W, int H, int pad) {
	return 0;
}

int removePadding(int* img, int W, int H, int pad) {
	return 0;
}

int kernel(int* neighbors, int nth) {
	int neighborsN = neighborsNumber(neighbors);
	if(nth % 2)
		return (neighborsN >= 3 && neighborsN <= 6) && (sequenceCounter(neighbors) == 1) && firstPassDiagonal(neighbors);
	else
		return (neighborsN >= 3 && neighborsN <= 6) && (sequenceCounter(neighbors) == 1) && secondPassDiagonal(neighbors);
		
}
*/
