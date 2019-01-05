#include <mpi.h>
#include <malloc.h>
#include <unistd.h>

#include "mpi_inst.h"
#include "utils.h"
#include "skeletonize.h"

#define getBlockIndex(block,W,i) (*I + block * W * i)

int myrank;

int mpi_init(int argc, char** argv) {
	MPI_Init(&argc, &argv);
	MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

	return 0;
}

static void mpi_ske_scatter(int **I, int block, int middle_block, int bottom_block, int W, int n_threads, int tag, MPI_Request *req) {
	int i;
	for(i=1; i < n_threads - 1; i++) { //we shall send 1 line above and 1 line below the block
		MPI_Isend( getBlockIndex(block,W,i) - W, middle_block * W, MPI_INT, i, tag, MPI_COMM_WORLD, req);
	}

	if(n_threads > 1) {
		// for the last process, we send a 1 pixel high slice (block-1) and the rest of the image
		MPI_Isend( getBlockIndex(block,W,i) - W, bottom_block * W, MPI_INT, i, tag, MPI_COMM_WORLD, req);
	}

}

static void mpi_ske_gather(int **I, int block, int H, int W, int n_threads, int gather_tag, MPI_Status* status) {
	// we receive blocks
	int i;
	for(i=1; i < n_threads - 1; i++) {
		MPI_Recv(getBlockIndex(block,W,i), block * W, MPI_INT, i, gather_tag, MPI_COMM_WORLD, status);
	}

	// last block to be received
	MPI_Recv(getBlockIndex(block,W,i), block * W + (H % n_threads) * W, MPI_INT, i, gather_tag, MPI_COMM_WORLD, status);
}

static void cleanup_padding(int *ch_image, int H, int W) {
	for(int u = 0; u < H-1; u++) {
		ch_image[u*W + 0] = 0;
		ch_image[u*W + (W-1)] = 0;
	}
	for(int u = 0; u < W-1; u++) {
		ch_image[0 + u] = 0;
		ch_image[(H-1)*W + u] = 0;
	}
}

int mpi_ske_start_centralized(int **I, int W, int H) {
	int n_threads = 1;
	int tag = 0;
	int source = 0;
	int i = 0;
	int* myimg, *ch_image, *aux;

	MPI_Status status;
	MPI_Request req;
	MPI_Comm_size(MPI_COMM_WORLD, &n_threads);

	int block = H/n_threads; //number of lines for each process
	int top_block = block + 1;
	int bottom_block = 1 + block + (H % n_threads);
	int middle_block = 1 + block + 1;

	/* scatter image from root into all other processes*/

	if(myrank == 0) {
		ch_image = (int*) memalign (32, H * W * sizeof(int)); 
		cleanup_padding(ch_image, H, W);

		mpi_ske_scatter(I,block,middle_block,bottom_block,W,n_threads,tag,&req);
	} else {
		// we receive blocks
		for(i=1; i < n_threads - 1; i++) {
			if(myrank == i) {
				myimg = (int*) memalign(0x20, middle_block * W * sizeof(int));
				ch_image = (int*) memalign (32, middle_block * W * sizeof(int)); 
				cleanup_padding(ch_image, middle_block, W);

				MPI_Recv(myimg , middle_block * W, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
			}
		}

		// last block to be received
		if(n_threads > 1 && myrank == n_threads - 1) {
			myimg = (int*) memalign(0x20, bottom_block * W * sizeof(int));
			ch_image = (int*) memalign (32, bottom_block * W * sizeof(int)); 
			cleanup_padding(ch_image, bottom_block, W);

			MPI_Recv(myimg , bottom_block * W, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
		}
	}

	// we work on them
	int cont0 = 1;
	int contOthers = 1;
	int iteration = 0;
	while(cont0 > 0) {
		if(n_threads == 1) { cont0 = skeletonize_matrixswap_dist(I, &ch_image, W, block, iteration); iteration++; continue; }
		if(myrank == 0) { // i'm with the top
			contOthers = skeletonize_matrixswap_dist(I, &ch_image, W, block+1, iteration);

			/* Gather */

			mpi_ske_gather(I,block,H,W,n_threads,tag,&status);

			/* Scatter */

			mpi_ske_scatter(I,block,middle_block,bottom_block,W,n_threads,tag,&req);
		}

		if(myrank != 0 && myrank != n_threads-1) { // middle
			contOthers = skeletonize_matrixswap_dist(&myimg, &ch_image, W, middle_block, iteration);
			MPI_Send( myimg + W, block * W, MPI_INT, source, tag, MPI_COMM_WORLD);
			MPI_Recv(myimg , middle_block * W, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
		}

		if(myrank == n_threads-1) { // i'm with the bottom part
			contOthers = skeletonize_matrixswap_dist(&myimg, &ch_image, W, bottom_block, iteration);
			MPI_Send( myimg + W, block * W + (H % n_threads) * W, MPI_INT, source, tag, MPI_COMM_WORLD);
			MPI_Recv(myimg , bottom_block * W, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
		}

		MPI_Allreduce(&contOthers, &cont0, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		//if(myrank == 0) { printf("cont0: %d\n", cont0); print_img(I,W,H); printf("\n"); }

		iteration++;
	}

	return iteration;
}

int mpi_ske_start(int **I, int W, int H) {
	int n_threads = 1;
	int tag = 2;
	int source = 0;
	int i = 0;
	int* myimg, *ch_image;
	int gather_tag = 1;

	MPI_Status status;
	MPI_Request req;
	MPI_Comm_size(MPI_COMM_WORLD, &n_threads);

	int block = H/n_threads; //number of lines for each process
	int top_block = block + 1;
	int bottom_block = 1 + block + (H % n_threads);
	int middle_block = 1 + block + 1;

	/* scatter image from root into all other processes*/

	if(myrank == 0) {
		ch_image = (int*) memalign (32, H * W * sizeof(int)); 
		cleanup_padding(ch_image, H, W);

		mpi_ske_scatter(I,block,middle_block,bottom_block,W,n_threads,tag, &req);
	} else {
		// we receive blocks
		for(i=1; i < n_threads - 1; i++) {
			if(myrank == i) {
				myimg = (int*) memalign(0x20, (middle_block * W) * sizeof(int));
				ch_image = (int*) memalign (32, middle_block * W * sizeof(int)); 
				cleanup_padding(ch_image, middle_block, W);

				MPI_Recv(myimg , middle_block * W, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
			}
		}

		// last block to be received
		if(n_threads > 1 && myrank == n_threads - 1) {
			myimg = (int*) memalign(0x20, bottom_block * W * sizeof(int));
			ch_image = (int*) memalign (32, bottom_block * W * sizeof(int)); 
			cleanup_padding(ch_image, bottom_block, W);

			MPI_Recv(myimg , bottom_block * W, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
		}
	}

	// we work on them
	int cont0 = 1;
	int contOthers = 1;
	int iteration = 0;
	while(cont0 > 0) {
		if(n_threads == 1) { cont0 = skeletonize_matrixswap_dist(I, &ch_image, W, block, iteration); iteration++; continue; } // jump
		if(myrank == 0) { // i'm with the top
			contOthers = skeletonize_matrixswap_dist(I, &ch_image, W, top_block, iteration);

			MPI_Isend(*I + (block-1)*W, W, MPI_INT, myrank+1, tag, MPI_COMM_WORLD, &req); 

			MPI_Recv(*I + block*W, W, MPI_INT, myrank+1, tag, MPI_COMM_WORLD, &status);

		}

		if(myrank != 0 && myrank != n_threads-1) { // middle
			contOthers = skeletonize_matrixswap_dist(&myimg, &ch_image, W, middle_block, iteration);

			MPI_Isend( myimg + W, W, MPI_INT, myrank-1, tag, MPI_COMM_WORLD, &req);
			MPI_Isend( &myimg[(middle_block-2)*W], W, MPI_INT, myrank+1, tag, MPI_COMM_WORLD, &req);

			MPI_Recv(myimg, W, MPI_INT, myrank-1, tag, MPI_COMM_WORLD, &status);
			MPI_Recv(&myimg[(middle_block - 1)*W] , W, MPI_INT, myrank+1, tag, MPI_COMM_WORLD, &status);

		}

		if(myrank == n_threads-1) { // i'm with the bottom part
			contOthers = skeletonize_matrixswap_dist(&myimg, &ch_image, W, bottom_block, iteration);

			MPI_Isend(myimg + W, W, MPI_INT, myrank-1, tag, MPI_COMM_WORLD,&req);

			MPI_Recv(myimg, W, MPI_INT, myrank-1, tag, MPI_COMM_WORLD, &status);
		}

		MPI_Allreduce(&contOthers, &cont0, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		iteration++;
	}

	/* Gather */

	if(n_threads > 1) {

		if(myrank == 0) { // I'm with the top 
			mpi_ske_gather(I,block,H,W,n_threads,gather_tag,&status);
		}

		if(myrank != 0 && myrank != n_threads-1) { // middle
			MPI_Ssend( myimg + W, block * W, MPI_INT, source, gather_tag, MPI_COMM_WORLD);
		}

		if(myrank == n_threads-1) { // I'm with the bottom part
			MPI_Ssend( myimg + W, block * W + (H % n_threads) * W, MPI_INT, source, gather_tag, MPI_COMM_WORLD);
		}

	}

	return iteration;
}

int mpi_ske_start_comm(int **I, int W, int H, int num_it) {
	int n_threads = 1;
	int tag = 2;
	int source = 0;
	int i = 0;
	int* myimg, *ch_image;
	int gather_tag = 1;

	MPI_Status status;
	MPI_Request req;
	MPI_Comm_size(MPI_COMM_WORLD, &n_threads);

	int block = H/n_threads; //number of lines for each process
	int top_block = block + 1;
	int bottom_block = 1 + block + (H % n_threads);
	int middle_block = 1 + block + 1;

	/* scatter image from root into all other processes*/

	if(myrank == 0) {
		ch_image = (int*) memalign (32, H * W * sizeof(int)); 
		cleanup_padding(ch_image, H, W);

		mpi_ske_scatter(I,block,middle_block,bottom_block,W,n_threads,tag,&req);
	} else {
		// we receive blocks
		for(i=1; i < n_threads - 1; i++) {
			if(myrank == i) {
				myimg = (int*) memalign(0x20, (middle_block * W) * sizeof(int));
				ch_image = (int*) memalign (32, middle_block * W * sizeof(int)); 
				cleanup_padding(ch_image, middle_block, W);

				MPI_Recv(myimg , middle_block * W, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
			}
		}

		// last block to be received
		if(n_threads > 1 && myrank == n_threads - 1) {
			myimg = (int*) memalign(0x20, bottom_block * W * sizeof(int));
			ch_image = (int*) memalign (32, bottom_block * W * sizeof(int)); 
			cleanup_padding(ch_image, bottom_block, W);

			MPI_Recv(myimg , bottom_block * W, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
		}
	}

	// we work on them
	int cont0 = 1;
	int contOthers = 1;
	int iteration = 0;
	while(cont0 > 0) {
		if(n_threads == 1) { cont0 = skeletonize_matrixswap_dist(I, &ch_image, W, block, iteration); iteration++; continue; } // jump
		if(myrank == 0) { // i'm with the top
			//contOthers = skeletonize_matrixswap_dist(I, &ch_image, W, top_block, iteration);
			num_it = num_it - 1;
			contOthers = num_it;

			if(contOthers == 0) { (*I)[(block-1)*W] = 1; }

			MPI_Isend(*I + (block-1)*W, W, MPI_INT, myrank+1, tag, MPI_COMM_WORLD, &req); 

			(*I)[(block-1)*W] = 0;

			MPI_Recv(*I + block*W, W, MPI_INT, myrank+1, tag, MPI_COMM_WORLD, &status);

			if((*I)[block*W] == 1) { (*I)[block*W] = 0; }
		}

		if(myrank != 0 && myrank != n_threads-1) { // middle
			//contOthers = skeletonize_matrixswap_dist(&myimg, &ch_image, W, middle_block, iteration);
			num_it = num_it - 1;
			contOthers = num_it;

			if(contOthers == 0) {	
				myimg[W] = 1;
				myimg[(middle_block-2)*W] = 1;
			} 

			MPI_Isend( myimg + W, W, MPI_INT, myrank-1, tag, MPI_COMM_WORLD, &req);
			MPI_Isend( &myimg[(middle_block-2)*W], W, MPI_INT, myrank+1, tag, MPI_COMM_WORLD, &req);

			myimg[W] = 0;
			myimg[(middle_block-2)*W] = 0;

			MPI_Recv(myimg, W, MPI_INT, myrank-1, tag, MPI_COMM_WORLD, &status);
			MPI_Recv(&myimg[(middle_block - 1)*W] , W, MPI_INT, myrank+1, tag, MPI_COMM_WORLD, &status);

			if(myimg[0] == 1) { myimg[0] = 0; }

			if(myimg[(block+1)*W] == 1) { myimg[(block+1)*W] = 0; }
		}

		if(myrank == n_threads-1) { // i'm with the bottom part
			//contOthers = skeletonize_matrixswap_dist(&myimg, &ch_image, W, bottom_block, iteration);
			num_it = num_it - 1;
			contOthers = num_it;

			if(contOthers == 0) { myimg[W] = 1; } 

			MPI_Isend(myimg + W, W, MPI_INT, myrank-1, tag, MPI_COMM_WORLD,&req);

			myimg[W] = 0;

			MPI_Recv(myimg, W, MPI_INT, myrank-1, tag, MPI_COMM_WORLD, &status);

			if(myimg[0] == 1) { myimg[0] = 0; }
		}

		MPI_Allreduce(&contOthers, &cont0, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		iteration++;
	}

	/* Gather */

	if(n_threads > 1) {

		if(myrank == 0) { // I'm with the top 
			mpi_ske_gather(I,block,H,W,n_threads,gather_tag,&status);
		}

		if(myrank != 0 && myrank != n_threads-1) { // middle
			MPI_Ssend( myimg + W, block * W, MPI_INT, source, gather_tag, MPI_COMM_WORLD);
		}

		if(myrank == n_threads-1) { // I'm with the bottom part
			MPI_Ssend( myimg + W, block * W + (H % n_threads) * W, MPI_INT, source, gather_tag, MPI_COMM_WORLD);
		}

	}

	return iteration;
}

int mpi_finalize() {
	//only the master thread should call this. In the event of using processes everyone should call it
	return MPI_Finalize();
}
