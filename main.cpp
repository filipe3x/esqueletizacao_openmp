#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <papi.h>
#include <omp.h>

#include "papi_inst.h"

#include "skeletonize_v3.h"
#include "ppm.h"
#include "utils.h"

static int verify_command_line (int argc, char *argv[], char *infile, char *outfile, int *fcode, int *f_width, int *num_threads);
static int read_inp_image (char *infile, image *img);
static int write_out_image (char *outfile, image res_img);
static void print_usage (char *msg);

int main (int argc, char *argv[]) {
  image img;
  int f_width, i, fcode, num_threads = 1;
  char infile[256], outfile[256];


  if (!verify_command_line (argc, argv, infile, outfile, &fcode, &f_width, &num_threads)) {
	return 0;
  }
    //printf("here is line %d\n", __LINE__); //for debugging


  // Read input file
  if (!read_inp_image (infile, &img)) return 0;

  // Initialize Papi and its events
  //int nbr_papi_runs = papi_init ();
  int nbr_papi_runs = 1; //for testing we just run one time
  if (!nbr_papi_runs) {
    fprintf (stderr, "Error initializing PAPI and corresponding events!\n");
    return 0;
  }

  for (i = 0; i < nbr_papi_runs; i++) {

     long long PAPI_start, PAPI_stop; 
     double stop, start = omp_get_wtime();
     PAPI_start = PAPI_get_real_usec();
     int it = 0;

     if(img->width < 55) print_img (img->buf, img->height, img->width); //print original image

     printf("Image size %d x %d \n", img->height, img->width); 
     printf("changing matrix allocated of size: %ld Kbytes\n", img->height * img->width * sizeof(int)/1024);

 //   papi_start_event (i);

     switch (fcode) {
       case 0:
	  //omp_set_num_threads(num_threads);
	  it = skeletonize(img->buf, img->width, img->height);
	  break;
       case 1:
	  it = skeletonize_serial(img->buf, img->width, img->height);
	  break;
       default:
	    print_usage ((char *)"Unknown function code!");
	    return 0;
	    break;
     }	 

//     papi_stop_event (i);

     if(img->width < 55) print_img (img->buf, img->height, img->width); //print result

     stop = omp_get_wtime();
     PAPI_stop = PAPI_get_real_usec();
     printf("papi Time in microseconds: %lld\n", PAPI_stop - PAPI_start);
     printf("omp Time in microseconds: %f\n",(stop - start)*1000000 );
     printf("total iterations: %d\n", it);
     printf("time per iteration in us: %.2f\n",  ( (double) (stop-start) * 1000000) /  (double) it );

  }
  
//  papi_print ();
  
//  papi_finalize ();


  if (!write_out_image (outfile, img)) return 0;

  free_img (img);

}

static int verify_command_line (int argc, char *argv[], char *infile, char *outfile, int *fcode, int *f_width, int *num_threads) {
	if (argc<4) {
		print_usage ((char *)"At least 3 arguments are required!");
		return 0;
	}
	strcpy (infile, argv[1]);
	strcpy (outfile, argv[2]);
	*fcode = atoi (argv[3]);
	switch (*fcode) {
	  case 0:
	    if (argc<5) {
		print_usage ((char *)"skeletonize parallel requires number of threads!");
		return 0;
	    }
	    *num_threads = atoi (argv[4]);
	    printf("running skeletonize parallel with %d threads\n", *num_threads);
	    break;
	  case 1:
	    printf("running skeletonize serial\n"); 
	    break;
	  default:
	    print_usage ((char *)"Unknown function code!");
	    return 0;
	    break;
	}
	return 1;
}

static void print_usage (char *msg) {
	fprintf (stderr, "Command Line Error! %s\n", msg);
	fprintf (stderr, "Usage:\tskeletonize <input filename> <output filename> <function code> [Specific Parameters]\n\n");
	fprintf (stderr, "\t<function code> = 0 : skeletonize_Parallel [Specific Parameters] = <n_threads>\n");
	fprintf (stderr, "\t<function code> = 1 : skeletonize_Serial [Specific Parameters] = NULL\n");
	fprintf (stderr, "\n");
}

static int read_inp_image (char *infile, image *img) {
	FILE *f;
	
	// open input file
	// NOTE: it must be open as binary, otherwise reading might stop at \0 chars
	f = fopen (infile, "rb");
	if (!f) {
		fprintf (stderr, "Error opening input file: %s\n", infile);
		return 0;
	}

	// read file
	*img = read_ppm(f);
	if (!(*img)) {
		fclose (f);
		fprintf (stderr, "Error reading image file: %s\n", infile);
		return 0;
	}
	fclose (f);
	return 1;
}

static int write_out_image (char *outfile, image res_img) {
	FILE *f;
	
	// open output file
	// NOTE: it must be open as binary, otherwise errors arise with the CR+LF and the LF convention of Win/Unix
	f = fopen (outfile, "wb");
	if (!f) {
		fprintf (stderr, "Error opening output file: %s\n", outfile);
		return 0;
	}
	output_ppm(f, res_img);
	fclose (f);

	return 1;
}
