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

static void mpi_ske_scatter(int *I, int block, int middle_block, int bottom_block, int W, int n_threads, int tag) {
	int i;
	for(i=1; i < n_threads - 1; i++) { //we shall send 1 line above and 1 line below the block
		MPI_Send( getBlockIndex(block,W,i) - W, middle_block * W, MPI_INT, i, tag, MPI_COMM_WORLD);
	}

	if(n_threads > 1) {
		// for the last process, we send a 1 pixel high slice (block-1) and the rest of the image
		MPI_Send( getBlockIndex(block,W,i) - W, bottom_block * W, MPI_INT, i, tag, MPI_COMM_WORLD);
	}

}

static void mpi_ske_gather(int *I, int block, int H, int W, int n_threads, MPI_Status* status) {
	// we receive blocks
	int i;
	for(i=1; i < n_threads - 1; i++) {
		MPI_Recv(getBlockIndex(block,W,i), block * W, MPI_INT, i, MPI_ANY_TAG, MPI_COMM_WORLD, status);
	}

	// last block to be received
	MPI_Recv(getBlockIndex(block,W,i), block * W + (H % n_threads) * W, MPI_INT, i, MPI_ANY_TAG, MPI_COMM_WORLD, status);
}

static int cleanup_padding(int *ch_image, int H, int W) {
	for(int u = 0; u < H-1; u++) {
		ch_image[u*W + 0] = 0;
		ch_image[u*W + (W-1)] = 0;
	}
	for(int u = 0; u < W-1; u++) {
		ch_image[0 + u] = 0;
		ch_image[(H-1)*W + u] = 0;
	}
}

int mpi_start_backup(int *I, int W, int H) {
	int n_threads = 1;
	int tag = 0;
	int source = 0;
	int i = 0;
	int* myimg, *ch_image, *aux;

	MPI_Status status;
	MPI_Comm_size(MPI_COMM_WORLD, &n_threads);

	int block = H/n_threads; //number of lines for each process
	int top_block = block + 1;
	int bottom_block = 1 + block + (H % n_threads);
	int middle_block = 1 + block + 1;

	/* scatter image from root into all other processes*/

	if(myrank == 0) {
		ch_image = (int*) memalign (32, H * W * sizeof(int)); 
		cleanup_padding(ch_image, H, W);

		mpi_ske_scatter(I,block,middle_block,bottom_block,W,n_threads,tag);
	} else {
		// we receive blocks
		for(i=1; i < n_threads - 1; i++) {
			if(myrank == i) {
				myimg = (int*) memalign(0x20, middle_block * W * sizeof(int));
				ch_image = (int*) memalign (32, middle_block * W * sizeof(int)); 
				cleanup_padding(ch_image, middle_block, W);

				MPI_Recv(myimg , middle_block * W, MPI_INT, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			}
		}

		// last block to be received
		if(n_threads > 1 && myrank == n_threads - 1) {
			myimg = (int*) memalign(0x20, bottom_block * W * sizeof(int));
			ch_image = (int*) memalign (32, bottom_block * W * sizeof(int)); 
			cleanup_padding(ch_image, bottom_block, W);

			MPI_Recv(myimg , bottom_block * W, MPI_INT, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		}
	}

	// we work on them
	int cont0 = 1;
	int contOthers = 1;
	int iteration = 0;
	while(cont0 > 0) {
		if(n_threads == 1) { cont0 = skeletonize_matrixswap_dist(&I, &ch_image, W, block, iteration); iteration++; continue; }
		if(myrank == 0) { // i'm with the top
			contOthers = skeletonize_matrixswap_dist(&I, &ch_image, W, block+1, iteration);

			/* Gather */

			mpi_ske_gather(I,block,H,W,n_threads,&status);

			/* Scatter */

			mpi_ske_scatter(I,block,middle_block,bottom_block,W,n_threads,tag);
		}

		if(myrank != 0 && myrank != n_threads-1) { // the virtue is in the middle
			contOthers = skeletonize_matrixswap_dist(&myimg, &ch_image, W, middle_block, iteration);
			MPI_Send( myimg + W, block * W, MPI_INT, source, tag, MPI_COMM_WORLD);
			MPI_Recv(myimg , middle_block * W, MPI_INT, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		}

		if(myrank == n_threads-1) { // i'm with the bottom part
			contOthers = skeletonize_matrixswap_dist(&myimg, &ch_image, W, bottom_block, iteration);
			MPI_Send( myimg + W, block * W + (H % n_threads) * W, MPI_INT, source, tag, MPI_COMM_WORLD);
			MPI_Recv(myimg , bottom_block * W, MPI_INT, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		}

		MPI_Allreduce(&contOthers, &cont0, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		//if(myrank == 0) { printf("cont0: %d\n", cont0); print_img(I,W,H); printf("\n"); }

		iteration++;
	}

	//MPI_Abort(MPI_COMM_WORLD, 999);

	return iteration;
}

int mpi_start(int *I, int W, int H) {
	int n_threads = 1;
	int tag = 0;
	int source = 0;
	int i = 0;
	int* myimg, *ch_image, *aux;

	MPI_Status status;
	MPI_Comm_size(MPI_COMM_WORLD, &n_threads);

	int block = H/n_threads; //number of lines for each process
	int top_block = block + 1;
	int bottom_block = 1 + block + (H % n_threads);
	int middle_block = 1 + block + 1;

	/* scatter image from root into all other processes*/

	if(myrank == 0) {
		ch_image = (int*) memalign (32, H * W * sizeof(int)); 
		cleanup_padding(ch_image, H, W);

		mpi_ske_scatter(I,block,middle_block,bottom_block,W,n_threads,tag);
	} else {
		// we receive blocks
		for(i=1; i < n_threads - 1; i++) {
			if(myrank == i) {
				myimg = (int*) memalign(0x20, (middle_block * W + 1) * sizeof(int));
				ch_image = (int*) memalign (32, middle_block * W * sizeof(int)); 
				cleanup_padding(ch_image, middle_block, W);

				MPI_Recv(myimg , middle_block * W, MPI_INT, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			}
		}

		// last block to be received
		if(n_threads > 1 && myrank == n_threads - 1) {
			myimg = (int*) memalign(0x20, bottom_block * W * sizeof(int));
			ch_image = (int*) memalign (32, bottom_block * W * sizeof(int)); 
			cleanup_padding(ch_image, bottom_block, W);

			MPI_Recv(myimg , bottom_block * W, MPI_INT, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		}
	}

	// we work on them
	int cont0 = 1;
	int contOthers = 1;
	int iteration = 0;
	int flag0 = 1;
	int flag1 = 1;
	while(contOthers > 0) {
		if(n_threads == 1) { cont0 = skeletonize_matrixswap_dist(&I, &ch_image, W, block, iteration); iteration++; continue; }
		if(myrank == 0) { // i'm with the top
			contOthers = skeletonize_matrixswap_dist(&I, &ch_image, W, top_block, iteration);

			//printf("it: %d, flag: %d\n", iteration, I[block*W]);

			if(contOthers == 0) { I[(block-1)*W] = 1; } // I didn't find any pixels, let's flag the neighbors to stop sending. My work is done

			if(flag1 == 1) MPI_Send(&I[(block-1)*W], W, MPI_INT, myrank+1, tag, MPI_COMM_WORLD); // We send our work if only they are still working as well. If they already finished, we stop bothering them. We must always signal our neighbors first before running away

			I[(block-1)*W] = 0;

			// Only if they are still working we wait and receive their messages
			if(flag1 == 1) MPI_Recv(&I[block*W], W, MPI_INT, myrank+1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			printf("contOthers: %d, rank: %d, it: %d, flag received: %d from rank %d\n", contOthers, myrank, iteration, I[block*W], myrank+1);

			// We received advice to stop bothering our neighbors. Let's respect that
			if(I[block*W] == 1) { flag1 = 0;  I[block*W] = 0; }

			printf("myrank: %d, myflag1: %d\n",myrank,flag1);

			// We keep doing our work nonetheless

		}


		// we assume if we finish our work, we don't send anything more, and our neighbors aren't expecting anything more from us. We just flag them, and we're out in our own way
		if(contOthers > 0 && myrank != 0 && myrank != n_threads-1) { // the virtue is in the middle
			contOthers = skeletonize_matrixswap_dist(&myimg, &ch_image, W, middle_block, iteration);

			if(contOthers == 0) { myimg[W] = 1; myimg[(middle_block-2)*W] = 1; } // I didn't find any pixels, let's flag the neighbors to stop sending. My work is done

			if(flag0 == 1) MPI_Send( myimg + W, W, MPI_INT, myrank-1, tag, MPI_COMM_WORLD); // We send our work if only they are still working as well. If they already finished, we stop bothering them. We must always signal our neighbors first before running away
			if(flag1 == 1) MPI_Send( &myimg[middle_block-2], W, MPI_INT, myrank+1, tag, MPI_COMM_WORLD);

			// Only if they are still working we wait and receive their messages
			if(flag0 == 1) MPI_Recv(myimg, W, MPI_INT, myrank-1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			printf("rank: %d, it: %d, flag0 received: %d from rank: %d\n", myrank, iteration, myimg[0],myrank-1);
			if(flag1 == 1) MPI_Recv(&myimg[(middle_block - 1)*W] , W, MPI_INT, myrank+1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			printf("rank: %d, it: %d, flag1 received: %d from rank: %d\n", myrank, iteration, myimg[(middle_block-1)*W], myrank+1);
			// We received advice to stop bothering our neighbors. Let's respect that
			if(myimg[0] == 1) { flag0 = 0; myimg[0] = 0; }
			if(myimg[(middle_block-1)*W] == 1) { flag1 = 0; myimg[(middle_block-1)*W] = 0; }

			// We keep doing our work nonetheless
		}

		if(contOthers > 0 && myrank == n_threads-1) { // i'm with the bottom part
			contOthers = skeletonize_matrixswap_dist(&myimg, &ch_image, W, bottom_block, iteration);

			if(contOthers == 0) { myimg[W] = 1; } // I didn't find any pixels, let's flag the neighbors to stop sending. My work is done

			if(flag0 == 1) MPI_Send(myimg + W, W, MPI_INT, myrank-1, tag, MPI_COMM_WORLD); // We send our work if only they are still working as well. If they already finished, we stop bothering them. We must always signal our neighbors first before running away

			myimg[W] = 0;

			// Only if they are still working we wait and receive their messages
			if(flag0 == 1) MPI_Recv(myimg, W, MPI_INT, myrank-1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

			printf("contOthers: %d, rank: %d, it: %d, flag received: %d from rank %d\n", contOthers, myrank, iteration, myimg[0],myrank-1);

			// We received advice to stop bothering our neighbors. Let's respect that
			if(myimg[0] == 1) { flag0 = 0; myimg[0] = 0; }

			// We keep doing our work nonetheless
		}

		//MPI_Allreduce(&contOthers, &contOthers, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);


		//if(myrank == 0) { printf("cont0: %d\n", cont0); print_img(I,W,H); printf("\n"); }

		iteration++;
	}
		if(myrank == 0) mpi_ske_gather(I,block,H,W,n_threads,&status);

		if(myrank != 0 && myrank != n_threads-1) { // the virtue is in the middle
			MPI_Send( myimg + W, block * W, MPI_INT, source, tag, MPI_COMM_WORLD);
		}

		if(myrank == n_threads-1) { // i'm with the bottom part
			MPI_Send( myimg + W, block * W + (H % n_threads) * W, MPI_INT, source, tag, MPI_COMM_WORLD);
		}

	//MPI_Abort(MPI_COMM_WORLD, 999);

	return iteration;
}

int mpi_finalize() {
	//only the master thread should call this. In the event of using processes everyone should call it
	return MPI_Finalize();
}
