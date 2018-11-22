#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

import sys
sys.path.append('/usr/lib/python3.6/site-packages')
from skimage.morphology import skeletonize as ske
from skimage import data
import matplotlib.pyplot as plt
from skimage.util import invert
import numpy as np

def createElipse(r):
        d = 2*r
        rx, ry = d/2, d/2
        x, y = np.indices((d, d))
        x, y = np.ogrid[-rx:d-rx, -ry:d-ry]
        mask = x*x + y*y <= rx*ry
        array = np.ones((d,d))
        array[mask] = 0
        return array.astype(np.uint8)

def im2double(im):
        min_val = np.min(np.ravel(im))
        max_val = np.max(np.ravel(im))
        out = (im.astype("float") - min_val) / (max_val - min_val)
        return out

def write_pgm(image, filename):
  """ Write grayscale image in PGM format to file.
  """
  height, width = image.shape
  maxval = image.max() ## should output 1 as max value
  with open(filename, 'wb') as f:
    f.write('P5\n{} {}\n{}\n'.format(width, height, 255))
    # not sure if next line works universally, but seems to work on my mac
    image.tofile(f)

def skeletonize(matrix, masksize):
  print "starting... matrix size: " + str(len(matrix[0])) + "; kernel size: " + str(masksize)
  kl = masksize
  ks = (kl-1)/2 ## kernels usually square with odd number of rows/columns
  img = matrix.copy();
  img = np.lib.pad(img,ks,'constant',constant_values=0) ## create padding of ks size with null pixels
  imx = len(img)
  imy = len(img[0])
  for k in to_infinity():
    print "it: " + str(k)
    list = []  ## initialize null list for saving pixels for deletion in this kth round
    for i in range(ks,imx-ks):
      for j in range(ks,imy-ks):
        mask = img[np.ix_([i-ks,i,i+ks],[j-ks,j,j+ks])] ## get mask/sub-matrix for i-j pixel of matrix
        if(img[i][j] == 1 and (border_elems(mask,1) == 0).sum() >= 1): ## we only apply masks to border elements
          #print "pixel: " + str(i) + "," + str(j)
          #if(i == 2 and j == 4): print mask
          if(kernel(mask,k % 2)): ## apply mask in this kth round
            #print "pixel: " + str(i) + "," + str(j)
            list.append((i,j)) ## mark pixel for deletion
    if(len(list) > 0): ## we have pixels to delete in this round
      for p in list: ## delete marked pixels
        #print "pixel removed: " + str(p)
        img.itemset(p,0) 
    elif(k % 2 == 1):
      break ## no more work to do
  return img[1:-1,1:-1] ## remove padding

def kernel(mask,nth):
        flag = 0
        if(neighborsNumber(mask) >= 2 and neighborsNumber(mask) <= 6):
                #print "number neig: " + str(neighborsNumber(mask))
                flag += 1
                #print "seq number:" + str(sequenceCounter(mask))
        if(sequenceCounter(mask) == 1):
                flag += 1
        if(nth == 0 and firstPassDiagonalNeighbors2(mask) == 0):
                #print "result mult first pass: " + str(firstPassDiagonalNeighbors(mask))
                flag += 1
        if(nth == 1 and secondPassDiagonalNeighbors2(mask) == 0):
                #print "result mult second pass: " + str(secondPassDiagonalNeighbors(mask))
                flag += 1
        if(flag == 3):
                return True
        else:
                return False

def neighborsNumber(mask):
        ones = (border_elems(mask,1) == 1).sum()
        return ones

def sequenceCounter(mask):
        seq = np.roll(border_elems(mask,1),-1)
        res = (seq[:-1] < seq[1:]).sum()
        if(seq[0] > seq[-1]): res = res + 1
        return res 

def firstPassDiagonalNeighbors(mask):
        return mask.item(5)^1 + mask.item(7)^1 + (mask.item(1)^1 * mask.item(3)^1)

def secondPassDiagonalNeighbors(mask):
        return mask.item(1)^1 + mask.item(3)^1 + (mask.item(5)^1 * mask.item(7)^1)

def firstPassDiagonalNeighbors2(mask):
        return mask.item(1) * mask.item(5) * mask.item(7) + mask.item(5) * mask.item(7) * mask.item(3) 

def secondPassDiagonalNeighbors2(mask):
        return mask.item(1) * mask.item(5) * mask.item(3) + mask.item(1) * mask.item(7) * mask.item(3)

def border_elems2(a, W): # Input array : a, Edgewidth : W
    n = a.shape[0]
    r = np.minimum(np.arange(n)[::-1], np.arange(n))
    return a[np.minimum(r[:,None],r)<W]

def border_elems(l, W): # Input array : l, Edgewidth : W
    r = list(l[0]) + list([i[-1] for i in l[1:-1]]) + list(reversed(l[-1])) + list(reversed([i[0] for i in l[1:-1]]))
    return np.asarray(r)

def to_infinity():
    index=0
    while 1:
        yield index
        index += 1

# Invert the horse image
X, Y = np.ogrid[0:9, 0:9]
image = (1./3 * (X - 4)**2 + (Y - 4)**2 < 3**2).astype(np.uint8)
image = invert(data.horse())
image = createElipse(767) ^ 1 ## real size will be (R+1)*2

# perform skeletonization
skeleton = image
#skeleton = ske(image)
#skeleton = skeletonize(image,3)

# save results
skeleton = np.lib.pad(skeleton,1,'constant',constant_values=0) ## create padding of ks size with null pixels
write_pgm(skeleton, "512Kcircle.pgm")

# display results
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(8, 4),
                         sharex=True, sharey=True)

ax = axes.ravel()

ax[0].imshow(image, cmap=plt.cm.gray)
ax[0].axis('off')
ax[0].set_title('original', fontsize=20)

ax[1].imshow(skeleton, cmap=plt.cm.gray)
ax[1].axis('off')
ax[1].set_title('skeleton', fontsize=20)

fig.tight_layout()
plt.show()

