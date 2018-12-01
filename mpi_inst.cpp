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

int mpi_start(int *I, int W, int H) {
	int n_threads = 1;
	int tag = 0;
	int source = 0;
	int i;
	int* myimg;

	MPI_Status status;
	MPI_Comm_size(MPI_COMM_WORLD, &n_threads);

	int block = H/n_threads; //number of lines for each process
	int top_block = block + 1;
	int bottom_block = 1 + block + (H % n_threads);
	int middle_block = 1 + block + 1;

	if(myrank == 0) {
		for(i=1; i < n_threads - 1; i++) { //we shall send 1 line above and 1 line below the block
			MPI_Ssend( getBlockIndex(block,W,i) - W, middle_block * W, MPI_INT, i, tag, MPI_COMM_WORLD);
		}
		printf("top block: \n");
		//print_img( I, W, block + 1 );
		printf("\n");

		if(n_threads > 1) {
			// for the last process, we send a 1 pixel high slice (block-1) and the rest of the image
			// print_img(W * (block - 1) * i + I, W, (block + 1) + (H % n_threads));
			MPI_Ssend( getBlockIndex(block,W,i) - W, bottom_block * W, MPI_INT, i, tag, MPI_COMM_WORLD);
			//printf("\nW: %d , H: %d\n", W,(block + 1) + (H % n_threads));
		}
	}


	// we receive the blocks
	for(i=1; i < n_threads - 1; i++) {
		if(myrank == i) {
			myimg = (int*) memalign(0x20, (block + 2) * W * sizeof(int)); // block+2 to account for the line above and below the block
			MPI_Recv(myimg , middle_block * W, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
			printf("middle block received: \n");
			//print_img( myimg, W, (block + 2) );
			printf("\n");
		}
	}

	// last block to be received
	if(n_threads > 1 && myrank == n_threads - 1) {
		myimg = (int*) memalign(0x20, bottom_block * W * sizeof(int));
		MPI_Recv(myimg , bottom_block * W, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
		printf("last block received: \n");
		//print_img( myimg, W, (block + 1) + (H % n_threads) );
		printf("\n");
	}


	//sleep(3);


	// we work on them
	int cont = 1;
	int passnr = 0;
	while(cont > 0 || passnr < 100) {
		if(n_threads == 1) { cont = skeletonize_matrixswap_dist(&I, W, block, passnr); printf("cont: %d, block: %d\n",cont,block); continue; }
		if(myrank == 0) { // i'm with the top
			printf("cont: %d, block: %d, pass: %d\n",cont,block,passnr);
			//print_img(I, W, H);
			printf("\n");
			cont = skeletonize_matrixswap_dist(&I, W, block+1, passnr);
			//printf("result rank0 cont: %d\n",cont);

			/* Gather and Scatter */

			// we receive blocks
			for(i=1; i < n_threads - 1; i++) {
				MPI_Recv(getBlockIndex(block,W,i), block * W, MPI_INT, i, tag, MPI_COMM_WORLD, &status);
			}
			// last block to be received
			MPI_Recv(getBlockIndex(block,W,i), block * W + (H % n_threads) * W, MPI_INT, i, tag, MPI_COMM_WORLD, &status);

			/* Scatter */

			for(i=1; i < n_threads - 1; i++) { //we shall send 1 line above and 1 line below the block
				MPI_Ssend( getBlockIndex(block,W,i) - W, middle_block * W, MPI_INT, i, tag, MPI_COMM_WORLD);
			}

			// for the last process, we send a 1 pixel high slice (block-1) and the rest of the image
			MPI_Ssend( getBlockIndex(block,W,i) - W, bottom_block * W, MPI_INT, i, tag, MPI_COMM_WORLD);
		}

		if(myrank != 0 && myrank != n_threads-1) { // the virtue is in the middle
			cont = skeletonize_matrixswap_dist(&myimg, W, middle_block, passnr);
			MPI_Ssend( myimg + W, block * W, MPI_INT, source, tag, MPI_COMM_WORLD);
			MPI_Recv(myimg , middle_block * W, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
		}

		if(myrank == n_threads-1) { // i'm with the bottom part
			cont = skeletonize_matrixswap_dist(&myimg, W, bottom_block, passnr);
			//printf("rank1:\n");
			//print_img(myimg, W, (block + 1) + (H % n_threads));
			//printf(" -- ^ rank1\n");
			MPI_Ssend( myimg + W, block * W + (H % n_threads) * W, MPI_INT, source, tag, MPI_COMM_WORLD);
			//print_img(myimg, W, (block + 1) + (H % n_threads));
			MPI_Recv(myimg , bottom_block * W, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
			//sleep(1);
		}

		passnr++;
		
	}

	//MPI_Abort(MPI_COMM_WORLD, 999);

	return 0;
}

int mpi_finalize() {
	//only the master thread should call this. In the event of using processes everyone should call it
	return MPI_Finalize();
}
