#include <mpi.h>

static int myrank;
static int n_threads;

int mpi_init(int t) {
	n_threads = t;
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
}

int mpi_start(int *I, int W, int H) {
	int block = H/n_threads; //number of lines for each of the hungry to treat
	int tag = 0;
	int source = 0;
	int i;
	int* myimg;

	if(myrank == 0) {
		for(i=1; i < n_threads - 1; i++) { //we shall send 1 line above and 1 line below the block
			MPI_Send( W * (block-1) * i + I, (block + 2) * W, MPI_INT, i, tag, MPI_COMM_WORLD);
		}

		// for the last process, we send a 1 pixel high slice (block-1) and the rest of the image
		MPI_Send( W * (block-1) * i + I, (block + 1) * W + (H % n_threads), MPI_INT, i, tag, MPI_COMM_WORLD);
	}

	// we receive the blocks
	for(i=1; i < n_threads - 1; i++) {
		if(myrank == i) {
			myimg = (int*) memalign(0x20, (block + 2) * W * sizeof(int)); // block+2 to account for the line above and below the block
			MPI_Recv(myimg , (block + 2) * W, MPI_INT, source, tag, MPI_COMM_WORLD);
		}
	}

	// last block to be received
	if(myrank == n_threads - 1) {
		myimg = (int*) memalign(0x20, ((block + 1) * W + H % n_threads) * sizeof(int)); // block+2 to account for the line above and below the block
		MPI_Recv(myimg , (block + 2) * W, MPI_INT, i, tag, MPI_COMM_WORLD);
	}

	// we work on them
	if(myrank == 0) { // i'm with the top
		skeletonize_matrixswap_serial(I, W, block+1);
	}

	if(myrank != 0 && myrank != n_threads-1) { // the virtue is in the middle
		skeletonize_matrixswap_serial(myimg, W, block+2);
	}

	if(myrank == n_threads-1) { // i'm with the bottom part
		skeletonize_matrixswap_serial(myimg, W, block+1);
	}
}

int mpi_finalize() {
	//only the master thread should call this. In the event of using processes everyone should call it
	return MPI_Finalize();
}
