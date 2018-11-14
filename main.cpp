#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <papi.h>

#include "papi_inst.h"

#include "skeletonize.h"
#include "ppm.h"
#include "utils.h"

static int verify_command_line (int argc, char *argv[], char *infile, char *outfile, int *fcode, int *f_width, int *num_threads);
static int read_inp_image (char *infile, image *img);
static int write_out_image (char *outfile, image res_img);
static void print_usage (char *msg);

//only gets used at the end to normalize the resultant matrix
void normalize (int *h, int W, int H) {
	int x, y;

	//normalize image
	for (y=0 ; y<H ; y++) { // for each row of I
		for (x=0 ; x<W ; x++) { // for each column of I
			if(h[y*W+x] > 255) h[y*W+x] = 255;
			if(h[y*W+x] < 0) h[y*W+x] = 0;
		}
	}

}

int main (int argc, char *argv[]) {
  image img, res_img, filter;
  int f_width, i, fcode, num_threads = 1;
  char infile[256], outfile[256];


  if (!verify_command_line (argc, argv, infile, outfile, &fcode, &f_width, &num_threads)) {
	return 0;
  }
    //printf("here is line %d\n", __LINE__); //for debugging


  // Read input file
  if (!read_inp_image (infile, &img)) return 0;

  // create output image
  res_img = new_img (img->width, img->height, BW);
  if (!res_img) {
	fprintf (stderr, "Error creating result image\n");
	return 0;
  } 

  // compute the filter for the desired width
  if (fcode==0)  {  // used only with convolve2D
    filter = new_img(f_width, f_width, BW);
    if (!filter) {
	fprintf (stderr, "Error creating filter\n");
	return 0;
    }
    if (!sharpen(filter->buf, f_width)) return 0;
    //print_gauss (filter->buf, f_width); //print the kernel used
  }

  // Initialize Papi and its events
  //int nbr_papi_runs = papi_init ();
  int nbr_papi_runs = 1; //for testing we just run one time
  if (!nbr_papi_runs) {
    fprintf (stderr, "Error initializing PAPI and corresponding events!\n");
    return 0;
  }

  for (i = 0; i < nbr_papi_runs; i++) {

//    long long PAPI_start, PAPI_stop; 
//    PAPI_start = PAPI_get_real_usec();

     print_gauss (img->buf, img->width); //print original image

 //   papi_start_event (i);

     switch (fcode) {
       case 0:     // convolve2D
	  skeletonize(res_img->buf, img->buf, img->width, img->height, filter->buf, f_width, num_threads);
	  break;
       case 1:
	  //convolve3x1 (res_img->buf, img->buf, img->width, img->height);
	  break;
       default:
	    print_usage ((char *)"Unknown function code!");
	    return 0;
	    break;
     }	 

//     papi_stop_event (i);

     print_gauss (res_img->buf, img->width); //print result

	//normalizing results
	//normalize (res_img->buf,  img->width, img->height);

//	PAPI_stop = PAPI_get_real_usec();
//	printf("Time in microseconds: %lld\n", PAPI_stop - PAPI_start);

  }
  
  papi_print ();
  
  papi_finalize ();


  if (!write_out_image (outfile, res_img)) return 0;

  free_img (filter); 
  free_img (img);
  free_img (res_img); 

  printf ("\nThat's all, folks\n");
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
	    if (argc<6) {
		print_usage ((char *)"skeletonize requires the filter_width parameter and number of threads!");
		return 0;
	    }
	    *f_width = atoi (argv[4]);
	    if (*(f_width)%2==0) {
		print_usage ((char *)"Filter width must be an odd number!");
		return 0;
	    }
		*num_threads = atoi (argv[5]);
	    break;
	  case 1:
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
	fprintf (stderr, "Usage:\tskeletonize <input filename> <output filename> <function code> [Specific Parameters] <filter_width> <n_threads>\n\n");
	fprintf (stderr, "\t<function code> = 0 : skeletonize [Specific Parameters] = <filter_width> <n_threads>\n");
	fprintf (stderr, "\t<function code> = 1 : convolve3x1 [Specific Parameters] = NULL\n");
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
}

// main function for adding noise to an image
// commented out!!!

/* void main (int argc, char *argv[]) {
	image img;
	FILE *f;
	int f_width;
	char infile[256], outfile[256];

	if (!verify_command_line (argc, argv, infile, outfile, &f_width)) {
		return;
	}

	// open input file
	// NOTE: it must be open as binary, otherwise reading might stop at \0 chars
	f = fopen (infile, "rb");
	if (!f) {
		fprintf (stderr, "Error opening input file: %s\n", infile);
		return;
	}

	img = read_ppm(f);
	if (!img) {
		fprintf (stderr, "Error reading image file:%s\n", infile);
		return;
	}
	fclose (f);

	add_noise (img, 100);
	// open output file
	// NOTE: it must be open as binary, otherwise errors arise with the CR+LF and the LF convention of Win/Unix
	f = fopen (outfile, "wb");
	if (!f) {
		fprintf (stderr, "Error opening output file: %s\n", outfile);
		return;
	}
	output_ppm(f, img);
	fclose (f);

	free_img (img); 


	printf ("\nThat's all, folks\n");
}
  */
