#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <papi.h>
#include <malloc.h>

#include "utils.h"
#include <iostream>
#include <cstdio>

using namespace std;

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
	int neighbors[10];									 // each pixel will have 8 neighbors
	int *ch_image = (int*) malloc (W * H * sizeof(int)); // which pixels we are going to change
	int *cont = (int*) malloc(2 * sizeof(int)); // for checking if we already finished

	/* Neighbors relative coordenates - P0 = total
						P0 P1 P2 P3 P4 P5 P6 P7 P8 P9 */
	int PX_index[10] = { 0, 0, 0, 1, 1, 1, 0,-1,-1,-1};
	int PY_index[10] = { 0, 0,-1,-1, 0, 1, 1, 1, 0,-1};

	int iterations = 0;		// total iterations count
	int i, j, k; 			// indexes
	int total, neigh_ant;	// total of neighbors and precious neighbor
	int ans;				// total transitions from 0 to 1
	int *im_read = I, *aux = NULL;

	cont[0] = 1;

	// cleaning the padding
	for(int u = 0; u < H-1; u++) {
		ch_image[u*W + 0] = BLACKPIXEL;
		ch_image[u*W + (W-1)] = BLACKPIXEL;
	}
	for(int u = 0; u < W-1; u++) {
		ch_image[0 + u] = BLACKPIXEL;
		ch_image[(H-1)*W + u] = BLACKPIXEL;
	}

	while((cont[0] > 0 || cont[1] > 0) || (iterations % 2 != 0)) {
		cont[0] = 0;
		cont[1] = 0;
		iterations = iterations + 1;
		//cout << "it=" << iterations << endl;

		for(i=1; i < H-1; i++) {
			for(j=1; j < W-1; j++) {
				/* error looking
				if((im_read[i*W + j] != 1) && (im_read[i*W + j] != 0))
					cout << im_read[i*W + j] << " em " << i << ", " << j << endl;*/
				// Only will do something if the point where it is at the moment is a 1;
				if(im_read[i*W + j] != BLACKPIXEL) {
					total = 0;
					ans = 0;
					neigh_ant = 1;	//redundante
					// get all neighbors
					for(k=1; k < 10; k++) {
						neighbors[k] = im_read[(i+PY_index[k])*W + (j+PX_index[k])];
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

					/*for(int u = 0; u < 10; u++)
						cout << neighbors[u] << ", ";
					cout << endl;*/

					if((total >= 2) && (total <= 6) && (ans == 1)) {
						if((iterations % 2 == 1) && (((1-neighbors[4]) + (1-neighbors[6]) + ((1-neighbors[2]) * (1-neighbors[8]))) >= 1)) {
							// we remove the pixel
							ch_image[i*W + j] = BLACKPIXEL;
							cont[0] += 1;
						} else {
							if((iterations % 2 == 0) && (((1-neighbors[2]) + (1-neighbors[8]) + ((1-neighbors[4]) * (1-neighbors[6]))) >= 1)) {
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
				// error looking
				if((ch_image[i*W + j] != 1) && (ch_image[i*W + j] != 0))
					cout << ch_image[i*W + j] << " em " << i << ", " << j << endl;
			}
		}

		// we delete pixels
		/* No we don't
		for(i=1; i < H-1; i++) {
			for(j=1; j < W-1; j++) {
				if(chan1to0[i*W+j] == 1) I[i*W+j] = 0;
			}
		}*/
		/*
		cout << "Imagem I\n";
		print_img(im_read, W, H);
		cout << "Imagem ch_image\n";
		print_img(ch_image, W, H);
		getchar();
		*/

		aux = im_read;
		im_read = ch_image;
		ch_image = aux;
	}
    return iterations;
}


// neighbors[k] = im_read[(i+PY_index[k])*W + (j+PX_index[k])];


int skeletonize (int *I, int W, int H) {
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

		// we delete pixels
		/* No we don't
		for(i=1; i < H-1; i++) {
			for(j=1; j < W-1; j++) {
				if(chan1to0[i*W+j] == 1) I[i*W+j] = 0;
			}
		}*/
		/*
		cout << "Imagem I\n";
		print_img(im_read, W, H);
		cout << "Imagem ch_image\n";
		print_img(ch_image, W, H);
		getchar();
		*/

		aux = im_read;
		im_read = ch_image;
		ch_image = aux;
	}
    return iterations;
}

/*
int skeletonize (int *I, int W, int H) {
	//int *neighbors = (int*) malloc(10 * sizeof(int)); // each pixel will have 8 neighbors
	int *ch_image = (int*) malloc (W * H * sizeof(int)); // which pixels we are going to change
	int *cont = (int*) malloc(2 * sizeof(int)); // for checking if we already finished

	// Neighbors relative coordenates - P0 = total
						P0 P1 P2 P3 P4 P5 P6 P7 P8 P9 //
	int PX_index[10] = { 0, 0, 0, 1, 1, 1, 0,-1,-1,-1};
	int PY_index[10] = { 0, 0,-1,-1, 0, 1, 1, 1, 0,-1};

	int iterations = 0;		// total iterations count
	//int i, j, k; 			// indexes
	int neigh_ant;	// total of neighbors and precious neighbor
	//int ans;				// total transitions from 0 to 1
	int *im_read = I, *aux = NULL;

	cont[0] = 1;

	// cleaning the padding
	for(int u = 0; u < H-1; u++) {
		ch_image[u*W + 0] = BLACKPIXEL;
		ch_image[u*W + (W-1)] = BLACKPIXEL;
	}
	for(int u = 0; u < W-1; u++) {
		ch_image[0 + u] = BLACKPIXEL;
		ch_image[(H-1)*W + u] = BLACKPIXEL;
	}

	while((cont[0] > 0 || cont[1] > 0) || (iterations % 2 != 0)) {
		cont[0] = 0;
		cont[1] = 0;
		iterations = iterations + 1;
		//cout << "it=" << iterations << endl;

		#pragma omp parallel for schedule(static) collapse(2)
		for(int i=1; i < H-1; i++) {
			for(int j=1; j < W-1; j++) {
				int neighbors[10], ans, total;
				// error looking
				if((im_read[i*W + j] != 1) && (im_read[i*W + j] != 0))
					cout << im_read[i*W + j] << " em " << i << ", " << j << endl;//
				// Only will do something if the point where it is at the moment is a 1;
				if(im_read[i*W + j] != BLACKPIXEL) {
					total = 0;
					ans = 0;
					neigh_ant = 1;	//redundante
					// get all neighbors
					for(int k=1; k < 10; k++) {
						neighbors[k] = im_read[(i+PY_index[k])*W + (j+PX_index[k])];
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

					//for(int u = 0; u < 10; u++)
						cout << neighbors[u] << ", ";
					cout << endl;//

					if((total >= 2) && (total <= 6) && (ans == 1)) {
						if((iterations % 2 == 1) && (((1-neighbors[4]) + (1-neighbors[6]) + ((1-neighbors[2]) * (1-neighbors[8]))) >= 1)) {
							// we remove the pixel
							ch_image[i*W + j] = BLACKPIXEL;
							cont[0] += 1;
						} else {
							if((iterations % 2 == 0) && (((1-neighbors[2]) + (1-neighbors[8]) + ((1-neighbors[4]) * (1-neighbors[6]))) >= 1)) {
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
				// error looking
				if((ch_image[i*W + j] != 1) && (ch_image[i*W + j] != 0))
					cout << ch_image[i*W + j] << " em " << i << ", " << j << endl;
			}
		}

		// we delete pixels
		
		//
		cout << "Imagem I\n";
		print_img(im_read, W, H);
		cout << "Imagem ch_image\n";
		print_img(ch_image, W, H);
		getchar();
		//

		aux = im_read;
		im_read = ch_image;
		ch_image = aux;
	}
    return iterations;
}
*/


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
