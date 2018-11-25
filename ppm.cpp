#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>

#include "ppm.h"

#define PPMREADBUFLEN 256

image read_ppm(FILE *pf) {
        char buf[PPMREADBUFLEN], *t;
        image img;
		unsigned char *temp_buf;
        unsigned int w, h, d;
        int r;
		ImageType type;
 
        if (pf == NULL) return NULL;
        t = fgets(buf, PPMREADBUFLEN, pf);
        /* the code fails if the white space following "P6" is not '\n' */
        if (t == NULL) return NULL;

		if ( !strncmp(buf, "P6\n", 3)  ) 
			type = COLOR;
		else if ( !strncmp(buf, "P5\n", 3)  ) 
			type = BW;
		else 
			return NULL;

        do
        { /* Px formats can have # comments after first line */
           t = fgets(buf, PPMREADBUFLEN, pf);
           if ( t == NULL ) return NULL;
        } while ( strncmp(buf, "#", 1) == 0 );
        r = sscanf(buf, "%u %u", &w, &h);
        if ( r < 2 ) return NULL;
 
        r = fscanf(pf, "%u", &d);
        if ( (r < 1) || ( d != 255 ) ) return NULL;
        fseek(pf, 1, SEEK_CUR); /* skip one byte, should be whitespace */
 
        img = new_img (w,h,type);
		if (!img) return NULL;

		temp_buf = (unsigned char *) memalign (0x80,w*h*(type == BW ? 1 : 3 ) * sizeof(unsigned char));
		if (!temp_buf) return NULL;

		size_t rd = fread(temp_buf, (type == BW ? 1 : 3 ) * sizeof(unsigned char), w*h, pf);
		if ( rd < w*h )  {
			free_img(img);
			free (temp_buf);
			return NULL;
		}

		// convert from temp_buf (unsigned char) to img->buf (int)
		unsigned char *from = temp_buf;
		int *to = img->buf;
		for (rd=0 ; rd < w*h ; rd++) {
			if (type==BW) {
				*to = (int)(*from);
				to++;
				from++;
			} else {   // COLOR
				*to = (int)luminance((float)(*(from)), (float)(*(from+1)), (float)(*(from+2)));
				to++;
				from+=3;
			}
		}
		free (temp_buf);
        return img;
}

void output_ppm(FILE *fd, image img)
{
  unsigned int n, i;
  int *from;
  unsigned char data;

  fprintf(fd, "P5\n%d %d\n1\n", img->width, img->height);
  n = img->width * img->height;
  for (i=0, from = img->buf ; i<n ; i++, from++) {
	  data = (unsigned char)(*from);
	  fwrite((void *)(&data), sizeof(unsigned char), 1, fd);
  }
  fflush(fd);
}
