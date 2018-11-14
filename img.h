#ifndef __IMG_H
#define __IMG_H

typedef enum {
	BW,
	COLOR
} ImageType;

typedef struct {
    unsigned int width;
    unsigned int height;
	ImageType type;
    int *buf;
} image_t;

typedef image_t * image; 

image new_img (int w, int h, ImageType t);
void free_img (image img);
void add_noise (image img, int amount);

#endif