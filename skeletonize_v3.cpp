#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <papi.h>
#include <malloc.h>

#define BLACKPIXEL 0
#define WHITEPIXEL 1

/*

Call:
void skeletonize(img->buf, img->width, img->height) 

The strategy is to aglomerate data and therefore, computations. Ideally we should have each data chunck fitting L2 cache. Each chunck should be given to a thread/core.
We can check if memory calls are flooding to L3 with papi counters.

Each pixel = 4 bytes

*/

int skeletonize_serial_v3 (int *I, int W, int H) {
	int *neighbors = (int*) malloc(10*sizeof(int)); // each pixel will have 8 neighbors
	int *ch_image = (int*) memalign (0x20, W*H*sizeof(int)); // which pixels we are going to change
	int *cont = (int*) malloc(2*sizeof(int)); // for checking if we already finished

	/* Neighbors relative coordenates - P0 = total
						P0 P1 P2 P3 P4 P5 P6 P7 P8 P9 */
	int PX_index[10] = { 0, 0, 0, 1, 1, 1, 0,-1,-1,-1};
	int PY_index[10] = { 0, 0,-1,-1, 0, 1, 1, 1, 0,-1};

	int iterations = 0;		// total iterations count
	int i, j, k; 			// indexes
	int total, neigh_ant;	// total of neighbors and precious neighbor
	int ans;				// total transitions from 0 to 1

	cont[0] = 1;

	while(cont[0] > 0 || cont[1] > 0) {
		cont[0] = 0;
		cont[1] = 0;
		iterations++;
		printf("it=%d\n",iterations);

		for(i=1; i < H-1; i++) {
			for(j=1; j < W-1; j++) {
				// Only will do something if the point where it is at the moment is a 1;
				if(I[i*W + j] == WHITEPIXEL) {
					total = 0;
					ans = 0;
					neigh_ant = 1;	//redundante
					// get all neighbors
					for(k=1; k < 10; k++) {
						neighbors[k] = I[(i+PY_index[k])*W + (j+PX_index[k])];
						if(k > 1) {
							total += neighbors[k];
							// previous IF could also be changed for a neighbors[0] = total-1; at end of cycle to take away the 1 from P1.
							if(neigh_ant - neighbors[k] == -1)
								ans++;
							// previous is the same as if(neigh_ant == 0 && neighbors[k] == 1)
						}
						neigh_ant = neighbors[k];
					}
					neighbors[0] = total;	// redundante for this part but might be important to paralelism
					if(neighbors[9] - neighbors[2] == -1)	// same as neigh_ant - neighbors[2]
						ans++;

					if((total >= 2) && (total <= 6) && (ans == 1)) {
						if((it % 2 != 0) && (neighbors[4] + neighbors[6] + neighbors[2] * neighbors[8] == 1)) {
							// we remove the pixel
							ch_image[i*W + j] = BLACKPIXEL;
							cont[0] += 1;
						} else {
							if((it % 2 == 0) && (neighbors[2] + neighbors[8] + neighbors[4] * neighbors[6] == 1)) {
								// we remove the pixel
								ch_image[i*W + j] = BLACKPIXEL;
								cont[1] += 1;
							} else {
								// we don't take it away
								ch_image[i*W + j] = WHITEPIXEL;
							}
						}

					} else {
						// We don't take the pixel away
						ch_image[i*W + j] = WHITEPIXEL;
					}
				} else {
					// only copying - sort of
					ch_image[i*W + j] = BLACKPIXEL;
				}
			}
		}

		// we delete pixels
		/* No we don't
		for(i=1; i < H-1; i++) {
			for(j=1; j < W-1; j++) {
				if(chan1to0[i*W+j] == 1) I[i*W+j] = 0;
			}
		}*/
	}

    return iterations;
}


int skeletonize (int *I, int W, int H) { 
	int *neighbors = (int*) malloc(9*sizeof(int)); // each pixel will have 8 neighbors
	int *chan1to0 = (int*) memalign (0x20, W*H*sizeof(int)); // which pixels we are going to change
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


int skeletonize_shared_neighbors (int *I, int W, int H) { 
	//int *neighbors = (int*) malloc(9*sizeof(int)); // each pixel will have 8 neighbors
	int t = omp_get_num_threads();
	int *neighbors = (int*) malloc(t*9*sizeof(int)); // each pixel will have 8 neighbors
	int *chan1to0 = (int*) memalign (0x20, W*H*sizeof(int)); // which pixels we are going to change
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

		#pragma omp target map (tofrom : cont[:2])
		{
			#pragma omp parallel for private(i, j, k, ans, total) schedule(static) collapse(1)
			for(i=1; i < H-1; i++) {
				int t_num = omp_get_thread_num();
				for(j=1; j < W-1; j++) {
					total = 0;
					ans = 0;
					chan1to0[i*W+j] = 0;

					for(k=0; k < 8; k++) { // get all neighbors
						neighbors[t*t_num+k] = I[(i+X_index[k])*W + (j+Y_index[k])];
						total += neighbors[t*t_num+k];
					}
					neighbors[8] = total;

					for(k=0; k < 7; k++) { // count nr transitions from 0 to 1
						if(neighbors[t*t_num+k] == 0 && neighbors[t*t_num+k+1] == 1)
							ans += 1;
					}

					if(neighbors[t*t_num+7] == 0 && neighbors[t*t_num+0] == 1) ans += 1;

					if(it % 2 != 0 && I[i*W+j] == 1 && neighbors[t*t_num+8] >= 2 && neighbors[t*t_num+8] <= 6 && ans == 1 
						&& neighbors[t*t_num+0] * neighbors[t*t_num+4] * neighbors[t*t_num+6] == 0 
						&& neighbors[t*t_num+2] * neighbors[t*t_num+4] * neighbors[t*t_num+6] == 0) {
						
						chan1to0[i*W+j] = 1; // we mark pixel for deletion
						cont[0] = 1;

					}

					if(it % 2 == 0 && I[i*W+j] == 1 && neighbors[t*t_num+8] >= 2 && neighbors[t*t_num+8] <= 6 && ans == 1 
						&& neighbors[t*t_num+0] * neighbors[t*t_num+2] * neighbors[t*t_num+4] == 0 
						&& neighbors[t*t_num+0] * neighbors[t*t_num+2] * neighbors[t*t_num+6] == 0) {
						
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
