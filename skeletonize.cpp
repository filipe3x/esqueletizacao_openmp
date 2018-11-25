#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <papi.h>
#include <malloc.h>

#define BLACKPIXEL 0
#define WHITEPIXEL 1

#define getNeighbor(i,j,k) I[(i+X_index[k])*W + (j+Y_index[k])]

/*

Call:
void skeletonize(img->buf, img->width, img->height) 

Each pixel = one int (4 bytes)

*/

int skeletonize_naive_serial (int *I, int W, int H) {
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
		//printf("it=%d\n",it = it + 1);

		for(i=1; i < H-1; i++) {
			for(j=1; j < W-1; j++) {
				total = 0;
				ans = 0;
				chan1to0[i*W+j] = 0;

				for(k=0; k < 8; k++) { // get all neighbors
					neighbors[k] = getNeighbor(i,j,k);
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


/* 
	Skeletonize,
	Double Pass version: For each iteration, we do a FIRST PASS to mark pixels and a SECOND PASS to delete them
	An aditional matrix int* chan1to0 is created to mark the pixels
*/
int skeletonize_doublepass_par (int *I, int W, int H) { 
	int t = omp_get_max_threads();
	int *chan1to0 = (int*) calloc (W*H,sizeof(int)); // which pixels we are going to change

	int X_index[8] = {-1,-1,0,1,1,1,0,-1}; // neighbors relative coordinates
	int Y_index[8] = {0,1,1,1,0,-1,-1,-1};

	int it = 0; // total iterations count
	int i, j, k; // indexes
	int total; // total of neighbors
	int ans; // total transitions from 0 to 1

	int cont0;
	int cont1;

	//a pool of threads is created here
	#pragma omp parallel shared(cont0, cont1) reduction(+:it) 
	{

	cont0 = 1;
	cont1 = 1;

	// every thread will test this condition
	while(cont0 > 0 || cont1 > 0) {

		#pragma omp barrier // mandatory to all threads be here before we set cont0 and cont1

		cont0 = 0;
		cont1 = 0;
		int flag0 = 0;
		int flag1 = 0;
		it = it + 1; //every thread will initialize it = 0, and will increment its own private it var, after we exit we broadcast everyones it value

		#pragma omp for private(i, j, k, ans, total) schedule(static)
		for(i=0; i < H-1; i++) {
			for(j=1; j < W-1; j++) {

				if(I[i*W+j] == 0){ chan1to0[i*W+j] = 0; continue; } // jump

				total = 0;
				ans = 0;
				chan1to0[i*W+j] = 0;

				for(k=0; k < 8; k++) { // count nr total neighbors
					total += getNeighbor(i,j,k);
				}

				for(k=0; k < 7; k++) { // count nr transitions from 0 to 1
					if(getNeighbor(i,j,k) == 0 && getNeighbor(i,j,k+1) == 1)
						ans += 1;
				}

				if(getNeighbor(i,j,7) == 0 && getNeighbor(i,j,0) == 1) ans += 1;

				if(it % 2 == 0 && total >= 2 && total <= 6 && ans == 1 
					&& getNeighbor(i,j,0) * getNeighbor(i,j,2) * getNeighbor(i,j,4) == 0 
					&& getNeighbor(i,j,0) * getNeighbor(i,j,2) * getNeighbor(i,j,6) == 0) {
					
					chan1to0[i*W+j] = 1; // we mark pixel for deletion
					flag1 = 1;

				}

				if(it % 2 != 0 && total >= 2 && total <= 6 && ans == 1 
					&& getNeighbor(i,j,0) * getNeighbor(i,j,4) * getNeighbor(i,j,6) == 0 
					&& getNeighbor(i,j,2) * getNeighbor(i,j,4) * getNeighbor(i,j,6) == 0) {
					
					chan1to0[i*W+j] = 1; // we mark pixel for deletion
					flag0 = 1;
				}
			}
		}

		/* implicit barrier */

		if(flag0 == 1) cont0 = 1;
		if(flag1 == 1) cont1 = 1;
	
		// we delete pixels
		#pragma omp for private(i,j) schedule(static) 
		for(i=1; i < H-1; i++) {
			#pragma omp simd
			for(j=1; j < W-1; j++) {
				I[i*W+j] = I[i*W+j] - chan1to0[i*W+j];
				//if(chan1to0[i*W+j] == 1) I[i*W+j] = 0; 
			}
		}

	}// while
	}// parallel block


	return it/t;
}


int skeletonize_doublepass_serial (int *I, int W, int H) { 
	int *chan1to0 = (int*) calloc (W*H,sizeof(int)); // which pixels we are going to change

	int X_index[8] = {-1,-1,0,1,1,1,0,-1}; // neighbors relative coordinates
	int Y_index[8] = {0,1,1,1,0,-1,-1,-1};

	int it = 0; // total iterations count
	int i, j, k; // indexes
	int total; // total of neighbors
	int ans; // total transitions from 0 to 1

	int cont0 = 1;
	int cont1;

	while(cont0 > 0 || cont1 > 0) {

		cont0 = 0;
		cont1 = 0;
		it = it + 1; 

		for(i=0; i < H-1; i++) {
			for(j=1; j < W-1; j++) {

				if(I[i*W+j] == 0){ chan1to0[i*W+j] = 0; continue; }

				total = 0;
				ans = 0;
				chan1to0[i*W+j] = 0;

				for(k=0; k < 8; k++) { // get all neighbors
					total += getNeighbor(i,j,k);
				}

				for(k=0; k < 7; k++) { // count nr transitions from 0 to 1
					if(getNeighbor(i,j,k) == 0 && getNeighbor(i,j,k+1) == 1)
						ans += 1;
				}

				if(getNeighbor(i,j,7) == 0 && getNeighbor(i,j,0) == 1) ans += 1;

				if(it % 2 != 0 && I[i*W+j] == 1 && total >= 2 && total <= 6 && ans == 1 
					&& getNeighbor(i,j,0) * getNeighbor(i,j,4) * getNeighbor(i,j,6) == 0 
					&& getNeighbor(i,j,2) * getNeighbor(i,j,4) * getNeighbor(i,j,6) == 0) {
					
					chan1to0[i*W+j] = 1; // we mark pixel for deletion
					cont0 = 1;
				}

				if(it % 2 == 0 && I[i*W+j] == 1 && total >= 2 && total <= 6 && ans == 1 
					&& getNeighbor(i,j,0) * getNeighbor(i,j,2) * getNeighbor(i,j,4) == 0 
					&& getNeighbor(i,j,0) * getNeighbor(i,j,2) * getNeighbor(i,j,6) == 0) {
					
					chan1to0[i*W+j] = 1; // we mark pixel for deletion
					cont1 = 1;

				}
			}
		}


		// we delete pixels
		for(i=1; i < H-1; i++) {
			#pragma omp simd
			for(j=1; j < W-1; j++) {
				I[i*W+j] = I[i*W+j] - chan1to0[i*W+j];
				//if(chan1to0[i*W+j] == 1) I[i*W+j] = 0; 
			}
		}

	}//while

	return it;
}

int skeletonize_matrixswap_par (int *I, int W, int H) {
	int *ch_image = (int*) memalign (32, W * H * sizeof(int)); // which pixels we are going to change
	//int *cont = (int*) malloc(2 * sizeof(int)); // for checking if we already finished

	/* Neighbors relative coordenates - P0 = total
						P0 P1 P2 P3 P4 P5 P6 P7 P8 P9 */
	int PX_index[10] = { 0, 0, 0, 1, 1, 1, 0,-1,-1,-1};
	int PY_index[10] = { 0, 0,-1,-1, 0, 1, 1, 1, 0,-1};

	int iterations = 0;		// total iterations count
	//int i, j, k; 			// indexes
	//int neigh_ant;	// total of neighbors and precious neighbor
	//int ans;				// total transitions from 0 to 1
	int *im_read = I, *aux = NULL;

	int cont0 = 1, cont1 = 1;

	// cleaning the padding
	#pragma omp parallel for
	for(int u = 0; u < H-1; u++) {
		ch_image[u*W + 0] = BLACKPIXEL;
		ch_image[u*W + (W-1)] = BLACKPIXEL;
	}
	#pragma omp parallel for
	for(int u = 0; u < W-1; u++) {
		ch_image[0 + u] = BLACKPIXEL;
		ch_image[(H-1)*W + u] = BLACKPIXEL;
	}

	while((cont0 > 0 || cont1 > 0) || (iterations % 2 != 0)) {
		cont0 = 0;
		cont1 = 0;
		iterations = iterations + 1;
		//cout << "it=" << iterations << endl;

		#pragma omp parallel for schedule(dynamic) reduction(+:cont0) reduction(+:cont1)
		for(int i=1; i < H-1; i++) {
			for(int j=1; j < W-1; j++) {
				int ans, total, neigh_ant;
				/* error looking
				if((im_read[i*W + j] != 1) && (im_read[i*W + j] != 0))
					cout << im_read[i*W + j] << " em " << i << ", " << j << endl;*/
				// Only will do something if the point where it is at the moment is a 1;
				if(im_read[i*W + j] != BLACKPIXEL) {
					total = 0;
					ans = 0;
					neigh_ant = 1;	//redundante
					// get all neighbors
					for(int k=2; k < 10; k++) {
						total += im_read[(i+PY_index[k])*W + (j+PX_index[k])];
						// previous IF could also be changed for a neighbors[0] = total-1; at end of cycle to take away the 1 from P1.
						if(neigh_ant - im_read[(i+PY_index[k])*W + (j+PX_index[k])] == -1)
							ans++;
						// previous is the same as if(neigh_ant == 0 && neighbors[k] == 1)
						neigh_ant = im_read[(i+PY_index[k])*W + (j+PX_index[k])];
					}
					if(im_read[(i+PY_index[9])*W + (j+PX_index[9])] - im_read[(i+PY_index[2])*W + (j+PX_index[2])] == -1)	// same as neigh_ant - neighbors[2]
						ans++;

					/*for(int u = 0; u < 10; u++)
						cout << neighbors[u] << ", ";
					cout << endl;*/

					if((total >= 2) && (total <= 6) && (ans == 1)) {
						if((iterations % 2 == 1) && (((1-im_read[(i+PY_index[4])*W + (j+PX_index[4])]) + (1-im_read[(i+PY_index[6])*W + (j+PX_index[6])]) + ((1-im_read[(i+PY_index[2])*W + (j+PX_index[2])]) * (1-im_read[(i+PY_index[8])*W + (j+PX_index[8])]))) >= 1)) {
							// we remove the pixel
							ch_image[i*W + j] = BLACKPIXEL;
							cont0 += 1;
						} else {
							if((iterations % 2 == 0) && (((1-im_read[(i+PY_index[2])*W + (j+PX_index[2])]) + (1-im_read[(i+PY_index[8])*W + (j+PX_index[8])]) + ((1-im_read[(i+PY_index[4])*W + (j+PX_index[4])]) * (1-im_read[(i+PY_index[6])*W + (j+PX_index[6])]))) >= 1)) {
								// we remove the pixel
								ch_image[i*W + j] = BLACKPIXEL;
								cont1 += 1;
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
				// error looking
				//if((ch_image[i*W + j] != 1) && (ch_image[i*W + j] != 0))
				//	cout << ch_image[i*W + j] << " em " << i << ", " << j << endl;
			}
		}

		aux = im_read;
		im_read = ch_image;
		ch_image = aux;
	}
    return iterations;
}

int skeletonize_matrixswap_serial (int *I, int W, int H) {
	int *ch_image = (int*) memalign (32, W * H * sizeof(int)); // which pixels we are going to change
	//int *cont = (int*) malloc(2 * sizeof(int)); // for checking if we already finished

	/* Neighbors relative coordenates - P0 = total
		P0 P1 P2 P3 P4 P5 P6 P7 P8 P9 */
	int PX_index[10] = { 0, 0, 0, 1, 1, 1, 0,-1,-1,-1};
	int PY_index[10] = { 0, 0,-1,-1, 0, 1, 1, 1, 0,-1};

	int iterations = 0;		// total iterations count
	//int i, j, k; 			// indexes
	//int neigh_ant;	// total of neighbors and precious neighbor
	//int ans;				// total transitions from 0 to 1
	int *im_read = I, *aux = NULL;

	int cont0 = 1, cont1 = 1;

	// cleaning the padding
	for(int u = 0; u < H-1; u++) {
		ch_image[u*W + 0] = BLACKPIXEL;
		ch_image[u*W + (W-1)] = BLACKPIXEL;
	}
	for(int u = 0; u < W-1; u++) {
		ch_image[0 + u] = BLACKPIXEL;
		ch_image[(H-1)*W + u] = BLACKPIXEL;
	}

	while((cont0 > 0 || cont1 > 0) || (iterations % 2 != 0)) {
		cont0 = 0;
		cont1 = 0;
		iterations = iterations + 1;

		for(int i=1; i < H-1; i++) {
			for(int j=1; j < W-1; j++) {
				int ans, total, neigh_ant;
				if(im_read[i*W + j] != BLACKPIXEL) {
					total = 0;
					ans = 0;
					neigh_ant = 1;	//redundante
					// get all neighbors
					for(int k=2; k < 10; k++) {
						total += im_read[(i+PY_index[k])*W + (j+PX_index[k])];
						// previous IF could also be changed for a neighbors[0] = total-1; at end of cycle to take away the 1 from P1.
						if(neigh_ant - im_read[(i+PY_index[k])*W + (j+PX_index[k])] == -1)
							ans++;
						// previous is the same as if(neigh_ant == 0 && neighbors[k] == 1)
						neigh_ant = im_read[(i+PY_index[k])*W + (j+PX_index[k])];
					}
					if(im_read[(i+PY_index[9])*W + (j+PX_index[9])] - im_read[(i+PY_index[2])*W + (j+PX_index[2])] == -1)	// same as neigh_ant - neighbors[2]
						ans++;

					if((total >= 2) && (total <= 6) && (ans == 1)) {
						if((iterations % 2 == 1) && (((1-im_read[(i+PY_index[4])*W + (j+PX_index[4])]) + (1-im_read[(i+PY_index[6])*W + (j+PX_index[6])]) + ((1-im_read[(i+PY_index[2])*W + (j+PX_index[2])]) * (1-im_read[(i+PY_index[8])*W + (j+PX_index[8])]))) >= 1)) {
							// we remove the pixel
							ch_image[i*W + j] = BLACKPIXEL;
							cont0 += 1;
						} else {
							if((iterations % 2 == 0) && (((1-im_read[(i+PY_index[2])*W + (j+PX_index[2])]) + (1-im_read[(i+PY_index[8])*W + (j+PX_index[8])]) + ((1-im_read[(i+PY_index[4])*W + (j+PX_index[4])]) * (1-im_read[(i+PY_index[6])*W + (j+PX_index[6])]))) >= 1)) {
								// we remove the pixel
								ch_image[i*W + j] = BLACKPIXEL;
								cont1 += 1;
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

		aux = im_read;
		im_read = ch_image;
		ch_image = aux;
	}
    return iterations;
}

