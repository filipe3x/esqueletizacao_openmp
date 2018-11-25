#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>

#include "img.h"

image new_img (int w, int h, ImageType t) {
	image img;
	img = (image) malloc (sizeof(image_t));
    if ( img == NULL ) return NULL;
	img->buf = (int *) memalign (0x80,w * h * sizeof(int));
    if ( img->buf == NULL ) return NULL;
	img->height = h;
	img->width = w;
	img->type = t;
	return img;
}

void free_img (image img) {
	free (img->buf);
	free (img);
}

void add_noise (image img, int amount) {
	int x,y;
	int noise, value, *ptr;

	for (ptr = img->buf, y=0 ; y< img->height ; y++) {
		for (x=0 ; x<img->width ; x++, ptr++) {
			noise = (int)(-1.f + 2.f * ((float)rand())/ ( RAND_MAX + 1.0 ));  // random number in [-1 .. 1[
			noise *= amount;
			value = *ptr + noise;
			if (value <= 0) *ptr = 0;
			else if (value >= 255) *ptr = 255;
			else *ptr = value;
		}
	}
}
