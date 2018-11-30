#include <mpi.h>
#include <malloc.h>

#include "mpi_inst.h"
#include "skeletonize.h"

int mpi_init(int argc, char** argv) {
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

	return 0;
}

int mpi_start(int *I, int W, int H) {
	int n_threads = 1;
	int tag = 0;
	int source = 0;
	int i;
	int* myimg;
	MPI_Status status;
	MPI_Comm_size(MPI_COMM_WORLD, &n_threads);

	int block = H/n_threads; //number of lines for each of the hungry to treat

	if(myrank == 0) {
		for(i=1; i < n_threads - 1; i++) { //we shall send 1 line above and 1 line below the block
			MPI_Send( W * (block-1) * i + I, (block + 2) * W, MPI_INT, i, tag, MPI_COMM_WORLD);
		}

		// for the last process, we send a 1 pixel high slice (block-1) and the rest of the image
		MPI_Send( W * (block-1) * i + I, (block + 1) * W + (H % n_threads) * W, MPI_INT, i, tag, MPI_COMM_WORLD);
	}

	// we receive the blocks
	for(i=1; i < n_threads - 1; i++) {
		if(myrank == i) {
			myimg = (int*) memalign(0x20, (block + 2) * W * sizeof(int)); // block+2 to account for the line above and below the block
			MPI_Recv(myimg , (block + 2) * W, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
		}
	}

	// last block to be received
	if(myrank == n_threads - 1) {
		myimg = (int*) memalign(0x20, ((block + 1) * W + H % n_threads) * sizeof(int));
		MPI_Recv(myimg , (block + 1) * W + (H % n_threads) * W, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
	}

	// we work on them
	int cont = 1;
	while(cont > 0) {
		if(myrank == 0) { // i'm with the top
			cont = skeletonize_matrixswap_dist(I, W, block+1);
			// we receive the blocks
			for(i=1; i < n_threads - 1; i++) {
				if(myrank == i) {
					MPI_Recv(I + (block * i) * W, block * W, MPI_INT, i, tag, MPI_COMM_WORLD, &status);
				}
			}
			// last block to be received
			if(myrank == n_threads - 1) {
				MPI_Recv(I + block * (n_threads - 1) * W, block * W + (H % n_threads) * W, MPI_INT, n_threads-1, tag, MPI_COMM_WORLD, &status);
			}
			for(i=1; i < n_threads - 1; i++) { //we shall send 1 line above and 1 line below the block
				MPI_Send( W * (block-1) * i + I, (block + 2) * W, MPI_INT, i, tag, MPI_COMM_WORLD);
			}

			// for the last process, we send a 1 pixel high slice (block-1) and the rest of the image
			MPI_Send( W * (block-1) * i + I, (block + 1) * W + (H % n_threads) * W, MPI_INT, i, tag, MPI_COMM_WORLD);
		}

		if(myrank != 0 && myrank != n_threads-1) { // the virtue is in the middle
			cont = skeletonize_matrixswap_dist(myimg, W, block+2);
			MPI_Send( myimg + W, block * W, MPI_INT, source, tag, MPI_COMM_WORLD);
			MPI_Recv(myimg , (block + 2) * W, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
		}

		if(myrank == n_threads-1) { // i'm with the bottom part
			cont = skeletonize_matrixswap_dist(myimg, W, block+1 + H % n_threads);
			MPI_Send( myimg + W, block * W + (H % n_threads) * W, MPI_INT, source, tag, MPI_COMM_WORLD);
		MPI_Recv(myimg , (block + 1) * W + (H % n_threads) * W, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
		}
	}

	return 0;
}

int mpi_finalize() {
	//only the master thread should call this. In the event of using processes everyone should call it
	return MPI_Finalize();
}
