#include <mpi.h>
#include <malloc.h>
#include <unistd.h>

#include "mpi_inst.h"
#include "utils.h"
#include "skeletonize.h"

#define getBlockIndex(block,W,i) (I + block * W * i)

int myrank;

int mpi_init(int argc, char** argv) {
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

	return 0;
}

void mpi_ske_scatter(int *I, int block, int middle_block, int bottom_block, int W, int n_threads, int tag) {
	int i;
	for(i=1; i < n_threads - 1; i++) { //we shall send 1 line above and 1 line below the block
		MPI_Ssend( getBlockIndex(block,W,i) - W, middle_block * W, MPI_INT, i, tag, MPI_COMM_WORLD);
	}

	if(n_threads > 1) {
		// for the last process, we send a 1 pixel high slice (block-1) and the rest of the image
		MPI_Ssend( getBlockIndex(block,W,i) - W, bottom_block * W, MPI_INT, i, tag, MPI_COMM_WORLD);
	}

}

void mpi_ske_gather(int *I, int block, int H, int W, int n_threads, MPI_Status* status) {
	// we receive blocks
	int i;
	for(i=1; i < n_threads - 1; i++) {
		MPI_Recv(getBlockIndex(block,W,i), block * W, MPI_INT, i, MPI_ANY_TAG, MPI_COMM_WORLD, status);
	}

	// last block to be received
	MPI_Recv(getBlockIndex(block,W,i), block * W + (H % n_threads) * W, MPI_INT, i, MPI_ANY_TAG, MPI_COMM_WORLD, status);
}

int mpi_start(int *I, int W, int H) {
	int n_threads = 1;
	int tag = 0;
	int source;
	int i;
	int* myimg, *ch_image;

	MPI_Status status;
	MPI_Comm_size(MPI_COMM_WORLD, &n_threads);

	int block = H/n_threads; //number of lines for each process
	int top_block = block + 1;
	int bottom_block = 1 + block + (H % n_threads);
	int middle_block = 1 + block + 1;

	/* scatter image from root into all other processes*/

	if(myrank == 0) {
		ch_image = (int*) memalign (32, H * W * sizeof(int)); 
		mpi_ske_scatter(I,block,middle_block,bottom_block,W,n_threads,tag);
	} else {
		// we receive blocks
		for(i=1; i < n_threads - 1; i++) {
			if(myrank == i) {
				myimg = (int*) memalign(0x20, middle_block * W * sizeof(int));
				ch_image = (int*) memalign (32, middle_block * W * sizeof(int)); 
				MPI_Recv(myimg , middle_block * W, MPI_INT, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			}
		}

		// last block to be received
		if(n_threads > 1 && myrank == n_threads - 1) {
			myimg = (int*) memalign(0x20, bottom_block * W * sizeof(int));
			ch_image = (int*) memalign (32, bottom_block * W * sizeof(int)); 
			MPI_Recv(myimg , bottom_block * W, MPI_INT, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		}
	}

	// we work on them
	int cont0 = 1;
	int contOthers = 1;
	int iteration = 0;
	while(cont0 > 0) {
		if(n_threads == 1) { cont0 = skeletonize_matrixswap_dist(I, ch_image, W, block, iteration); iteration++; continue; }
		if(myrank == 0) { // i'm with the top
			contOthers = skeletonize_matrixswap_dist(I, ch_image, W, block+1, iteration);

			/* Gather */

			mpi_ske_gather(I,block,H,W,n_threads,&status);

			/* Scatter */

			mpi_ske_scatter(I,block,middle_block,bottom_block,W,n_threads,tag);
		}

		if(myrank != 0 && myrank != n_threads-1) { // the virtue is in the middle
			contOthers = skeletonize_matrixswap_dist(myimg, ch_image, W, middle_block, iteration);
			MPI_Ssend( myimg + W, block * W, MPI_INT, source, tag, MPI_COMM_WORLD);
			MPI_Recv(myimg , middle_block * W, MPI_INT, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		}

		if(myrank == n_threads-1) { // i'm with the bottom part
			contOthers = skeletonize_matrixswap_dist(myimg, ch_image, W, bottom_block, iteration);
			MPI_Ssend( myimg + W, block * W + (H % n_threads) * W, MPI_INT, source, tag, MPI_COMM_WORLD);
			MPI_Recv(myimg , bottom_block * W, MPI_INT, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		}


		MPI_Allreduce(&contOthers, &cont0, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		//if(myrank == 0) printf("cont0: %d\n", cont0);

		iteration++;
	}

	//MPI_Abort(MPI_COMM_WORLD, 999);

	return iteration;
}

int mpi_finalize() {
	//only the master thread should call this. In the event of using processes everyone should call it
	return MPI_Finalize();
}
