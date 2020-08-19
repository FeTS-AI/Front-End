#include <iostream>
#include <cstdio>
#include <string>
#include "lddmm_common.h"
#include "lddmm_data.h"

int usage()
{
  printf("macf_pyramid: Paul's MACF implementation\n");
  printf("Usage: \n");
  printf("  macf_pyramid DIM NLEVELS sigma basename extension\n");
  printf("Parameters: \n"); 
  printf("  DIM      :  Image dimensionality\n");        // 1
  printf("  NLEVELS  :  Number of multi-resolution levels\n");        // 1
  printf("  sigma    :  Gaussian smoothing factor\n");        // 1
  printf("  basename :  Name of warps (all letters before Warp/InverseWarp)\n");        // 1
  printf("  extension:  Filename extension for input and output (nii.gz)\n");
  return -1;
}

template <class TFloat, uint VDim>
int my_main(int argc, char *argv[])
{
  typedef LDDMMData<TFloat, VDim> LDDMM;
  typedef typename LDDMM::ImageType ImageType;
  typedef typename LDDMM::VectorImageType VectorImageType;
  typedef typename LDDMM::ImagePointer ImagePointer;
  typedef typename LDDMM::VectorImagePointer VectorImagePointer;

  // No parameter checking!
  if(argc < 6) return usage();

  // Basic info
  unsigned int nlev = (unsigned int) atoi(argv[2]);
  double sigma = atof(argv[3]);
  char *basename = argv[4];
  char *extension = argv[5];

  // Read the warps
  char fnwarp[1024], fninvw[1024];
  sprintf(fnwarp, "%sWarp.%s", basename, extension);
  sprintf(fninvw, "%sInverseWarp.%s", basename, extension);

  VectorImagePointer w, iw;
  LDDMM::vimg_read(fnwarp, w);
  LDDMM::vimg_read(fninvw, iw);

  // For each level, resample image
  for (uint i = 1; i <= nlev; i++)
    {
    VectorImagePointer q, s;
    uint factor = 1 << i;

    char fnout[1024];
    sprintf(fnout, "%s_mr%d_Warp.%s", basename, i, extension);
    LDDMM::vimg_gaussian_smooth(w, false, factor * sigma, q);
    LDDMM::vimg_shrink(q, factor, s);
    LDDMM::vimg_clip_to_zero_in_place(s, 1.0e-4);
    LDDMM::vimg_write(s, fnout);

    sprintf(fnout, "%s_mr%d_InverseWarp.%s", basename, i, extension);
    LDDMM::vimg_gaussian_smooth(iw, false, factor * sigma, q);
    LDDMM::vimg_shrink(q, factor, s);
    LDDMM::vimg_clip_to_zero_in_place(s, 1.0e-4);
    LDDMM::vimg_write(s, fnout);
    }

  return 0;
}

int main(int argc, char *argv[])
{
  // Read input arguments
  if(argc < 2)
    return usage();

  // Get the number of dimensions
  uint ndim = (uint) atoi(argv[1]);
  if(ndim == 2)
    my_main<myreal, 2>(argc, argv);
  else if(ndim == 3)
    my_main<myreal, 3>(argc, argv);
  else
    return usage();

}
