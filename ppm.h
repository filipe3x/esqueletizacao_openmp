#ifndef __PPM
#define __PPM

#include "img.h"

image read_ppm(FILE *pf);

void output_ppm(FILE *fd, image img);

#define luminance(R,G,B) (0.2126f*R+0.7152f*G+0.0722f*B)

#endif