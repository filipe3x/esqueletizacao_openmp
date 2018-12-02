#ifndef __SKELE
#define __SKELE

int skeletonize_doublepass_serial (int *I, int W, int H);
int skeletonize_doublepass_par (int *I, int W, int H);
int skeletonize_matrixswap_serial (int *I, int W, int H);
int skeletonize_matrixswap_par (int *I, int W, int H);
int skeletonize_matrixswap_dist (int **I, int **ch_image, int W, int H, int passnr);
//extern inline int skeletonize_matrixswap_dist (int **I, int W, int H, int passnr) __attribute__((always_inline));

#endif
