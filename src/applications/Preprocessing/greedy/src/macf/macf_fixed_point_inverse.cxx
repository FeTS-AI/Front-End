#include <iostream>
#include <cstdio>
#include <string>
#include "lddmm_common.h"
#include "lddmm_data.h"

int usage()
{
  printf("macf_fixed_point_inverse: Paul's MACF implementation\n");
  printf("Usage: \n");
  printf("  macf_fixed_point_inverse DIM inwarp.nii outwarp.nii\n");
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
  if(argc < 4) return usage();

  // Basic info
  char *fninwarp = argv[2];
  char *fnoutwarp = argv[3];

  // Read the warp
  VectorImagePointer u, v, vn;
  LDDMM::vimg_read(fninwarp, u);

  // Allocate the inverse
  LDDMM::alloc_vimg(v, u);
  LDDMM::alloc_vimg(vn, u);

  // Perform fixed point algorithm
  for (uint i = 1; i <= 100; i++)
    {
    // Apply iteration
    LDDMM::interp_vimg(u, v, 1.0, true, vn);
    LDDMM::vimg_scale_in_place(vn, -1.0);

    // Measure change
    LDDMM::vimg_subtract_in_place(v, vn);
    myreal normdiff = LDDMM::vimg_euclidean_norm_sq(v);

    // Report 
    printf("Iteration %03d, Change %f\n", i, normdiff);
    LDDMM::vimg_scale(vn, 1.0, v);
    }

  LDDMM::vimg_write(v, fnoutwarp);

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
