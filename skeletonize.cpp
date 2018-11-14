#include <stdio.h>
#include <omp.h>
#include <papi.h>

/*

void skeletonize(res_img->buf, img->buf, img->width, img->height, filter->buf, f_width, num_threads)

The strategy is to aglomerate data and therefore, computations. Ideally we should have each data chunck fitting L2 cache. Each chunck should be given to a thread/core.
We can check if memory calls are flooding to L3 with papi counters.

*/
void skeletonize(int *h, int *I, int W, int H, int *f, int U, int num_threads) {
	int x, y, i, j;
	int halfU;
	int sumW;

	long long PAPI_start, PAPI_stop; 
	double stop, start = omp_get_wtime();
	PAPI_start = PAPI_get_real_usec();

	halfU = U/2;

	omp_set_num_threads(num_threads);
	// each thread will have its own value of x, sumW, i, j
	#pragma omp parallel for private(x, sumW, i, j) schedule(static)
	for (y=0 ; y<H ; y++) { // for each row of I
		for (x=0 ; x<W ; x++) { // for each column of I

			// compute h[y][x]
			sumW = h[y*W+x] = 0;
			for (i=-halfU ; i<=halfU ; i++) {
				// verify I horizontal bounds
				if (x+i<0 || x+i>=W) continue;

				for (j=-halfU ; j<=halfU ; j++) {
					// verify I vertical bounds
					if (y+j<0 || y+j>=H) continue;

					h[y*W+x] += (f[(j+halfU)*U+i+halfU] * I[(y+j)*W + (x+i)]);
					sumW += f[(j+halfU)*U+i+halfU];
				}
			}
			h[y*W+x] /= (sumW ? sumW : 1);
		} // y loop
	} // x loop

	stop = omp_get_wtime();
	PAPI_stop = PAPI_get_real_usec();
	printf("papi Time in microseconds: %lld\n", PAPI_stop - PAPI_start);
	printf("omp Time in microseconds: %f\n",(stop - start)*1000000 );
}

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

