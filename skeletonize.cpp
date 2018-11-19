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

Each pixel = 4 bytes

*/

int skeletonize_serial (int *I, int W, int H) {
	int *neighbors = (int*) malloc(9*sizeof(int)); // each pixel will have 8 neighbors
	int *chan1to0 = (int*) malloc(W*H*sizeof(int)); // which pixels we are going to change
	int *cont = (int*) malloc(2*sizeof(int)); // for checking if we already finished

	int X_index[8] = {-1,-1,0,1,1,1,0,-1}; // neighbors relative coordinates
	int Y_index[8] = {0,1,1,1,0,-1,-1,-1};

	int it = 0; // total iterations count
	int i, j, k; // indexes
	int total; // total of neighbors
	int ans; // total transitions from 0 to 1

	cont[0] = 1;

	while(cont[0] > 0 || cont[1] > 0) {
		cont[0] = 0;
		cont[1] = 0;
		it = it + 1;
		//printf("it=%d\n",it = it + 1);

		for(i=1; i < H-1; i++) {
			for(j=1; j < W-1; j++) {
				total = 0;
				ans = 0;
				chan1to0[i*W+j] = 0;

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

				if(it % 2 != 0 && I[i*W+j] == 1 && neighbors[8] >= 2 && neighbors[8] <= 6 && ans == 1 
					&& neighbors[0] * neighbors[4] * neighbors[6] == 0 
					&& neighbors[2] * neighbors[4] * neighbors[6] == 0) {
					
					chan1to0[i*W+j] = 1; // we mark pixel for deletion
					cont[0] = 1;

				}

				if(it % 2 == 0 && I[i*W+j] == 1 && neighbors[8] >= 2 && neighbors[8] <= 6 && ans == 1 
					&& neighbors[0] * neighbors[2] * neighbors[4] == 0 
					&& neighbors[0] * neighbors[2] * neighbors[6] == 0) {
					
					chan1to0[i*W+j] = 1; // we mark pixel for deletion
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

        return it;
}


int skeletonize (int *I, int W, int H) { 
	int *neighbors = (int*) malloc(9*sizeof(int)); // each pixel will have 8 neighbors
	int *chan1to0 = (int*) malloc(W*H*sizeof(int)); // which pixels we are going to change
	int *cont = (int*) malloc(2*sizeof(int)); // for checking if we already finished

	int X_index[8] = {-1,-1,0,1,1,1,0,-1}; // neighbors relative coordinates
	int Y_index[8] = {0,1,1,1,0,-1,-1,-1};

	int it = 0; // total iterations count
	int i, j, k; // indexes
	int total; // total of neighbors
	int ans; // total transitions from 0 to 1

	cont[0] = 1;

	while(cont[0] > 0 || cont[1] > 0) {
		cont[0] = 0;
		cont[1] = 0;
		it = it + 1;

		#pragma omp proc_bind(close)
		#pragma omp target map (to : neighbors[:9], X_index[:8], Y_index[:8]) map (tofrom : I[:W*H], chan1to0[:W*H], cont[:2])
		{
			#pragma omp parallel for private(i, j, k, ans, total) schedule(static) collapse(1)
			for(i=1; i < H-1; i++) {
				for(j=1; j < W-1; j++) {
					total = 0;
					ans = 0;
					chan1to0[i*W+j] = 0;

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

					if(it % 2 != 0 && I[i*W+j] == 1 && neighbors[8] >= 2 && neighbors[8] <= 6 && ans == 1 
						&& neighbors[0] * neighbors[4] * neighbors[6] == 0 
						&& neighbors[2] * neighbors[4] * neighbors[6] == 0) {
						
						chan1to0[i*W+j] = 1; // we mark pixel for deletion
						cont[0] = 1;

					}

					if(it % 2 == 0 && I[i*W+j] == 1 && neighbors[8] >= 2 && neighbors[8] <= 6 && ans == 1 
						&& neighbors[0] * neighbors[2] * neighbors[4] == 0 
						&& neighbors[0] * neighbors[2] * neighbors[6] == 0) {
						
						chan1to0[i*W+j] = 1; // we mark pixel for deletion
						cont[1] = 1;

					}
				}
			}

			// ** implicit barrier **

			// we delete pixels
			#pragma omp parallel for private(i, j) schedule(static) collapse(1)
			for(i=1; i < H-1; i++) {
				for(j=1; j < W-1; j++) {
					if(chan1to0[i*W+j] == 1) I[i*W+j] = 0;
				}
			}
		}
	}

	return it;
}

