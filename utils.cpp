// see math.h for explanation of the macro below (MS Visual Studio 2008)
#define _USE_MATH_DEFINES

#include <math.h>
#include <stdio.h>

int sharpen (int *f, int U) {
	int x,y, halfU;
	float sigma_sq, min_value, temp;

	if (U%2 == 0) {
		fprintf (stderr, "Filter width must be odd\n!");
		return 0;
	}

	halfU = U/2;

	for (y=-halfU ; y<=halfU ; y++) {
		for (x=-halfU ; x<=halfU ; x++, f++) {
			*f = 0;
			if(y == 0 || x == 0) *f = -1;
			if(y == 0 && x == 0) *f = U*2-1;

		}
	} 
	return 1;
}

int identity (int *f, int U) {
	int x,y, halfU;
	float sigma_sq, min_value, temp;

	if (U%2 == 0) {
		fprintf (stderr, "Filter width must be odd\n!");
		return 0;
	}

	halfU = U/2;

	for (y=-halfU ; y<=halfU ; y++) {
		for (x=-halfU ; x<=halfU ; x++, f++) {
			*f = 0;
			if(y == 0 && x == 0) *f = 1;

		}
	} 
	return 1;
}

int gauss (int *f, int U) {
	int x,y, halfU;
	float sigma_sq, min_value, temp;

	if (U%2 == 0) {
		fprintf (stderr, "Filter width must be odd\n!");
		return 0;
	}

	halfU = U/2;

	sigma_sq =((float)(U*U)-1.f)/12.f;

	// the [0][0] corner is the minimum value. Store it to convert all f elements to integers
	min_value = exp(-((float)(2*halfU*halfU))/(2.f*sigma_sq)) / (2.f * M_PI * sigma_sq);
	for (y=-halfU ; y<=halfU ; y++) {
		for (x=-halfU ; x<=halfU ; x++, f++) {
			temp = exp(-((float)(x*x+y*y))/(2.f*sigma_sq)) / (2.f * M_PI * sigma_sq);
			*f = (int)(temp/min_value);
		}
	}
	return 1;
}

void print_gauss (int *f, int U) {
	int j,i;

	for (j=0 ; j<U ; j++) {
		for (i=0 ; i<U ; i++) {
			fprintf (stdout, "%d\t", *f);
			f++;
		}
		fprintf (stdout, "\n");
	}
}

void print_img (int *f, int W, int H) {
	int j,i;

	for (j=0 ; j<H ; j++) {
		for (i=0 ; i<W ; i++) {
			fprintf (stdout, "%d ", f[W*j + i]);
		}
		fprintf (stdout, "\n");
	}
}


int addPadding(int* img, int W, int H, int pad) {
	return 0;
}

int removePadding(int* img, int W, int H, int pad) {
	return 0;
}
