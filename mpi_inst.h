
int mpi_init(int argc, char** argv);
int mpi_ske_start(int **I, int W, int H);
int mpi_ske_start_comm(int **I, int W, int H, int num_it);
int mpi_finalize();

extern int myrank;
