#include <iostream>
#include <cstdio>
#include "lddmm_common.h"
#include "lddmm_data.h"

/**
 * This program requires a directory to be in place. These are the files that will be accessed 
 * by the program:
 *
 * IDLIST.txt                          List of IDS (strings identifying subjects)
 * *_macf_InverseWarp.nii              Current estimate of warps from template to each image
 * fx_${ID}_mv_*_regWarp.nii           Warps from all images to image ID
 * fx_${ID}_mv_*_regInverseWarp.nii    Inverse warps from all images to image ID
 * fx_${ID}_mv_*_rank.nii              Weights associated with these warps
 * fx_*_mv_${ID}_regWarp.nii           Warps from image ID to all images
 * fx_*_mv_${ID}_regInverseWarp.nii    Inverse warps from image ID to all images
 * fx_*_mv_${ID}_rank.nii              Weights associated with these warps
 */

int usage()
{
  printf("macf_objective: Paul's MACF implementation\n");
  printf("Usage: \n");
  printf("  macf_objgrad DIM i\n");
  printf("Options: \n");
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

  typedef std::vector<ImagePointer> ImageArray;
  typedef std::vector<VectorImagePointer> VectorImageArray;

  // No parameter checking!
  if(argc < 4)
    return usage();

  // Basic info
  int N = atoi(argv[2]);
  int k = atoi(argv[3]);
  int iarg = 4;

  if(N <= 0 || i < 0 || i >= N) 
    return usage();

  if(argc < 4 + (N-1) * 5 + 2)
    return usage();

  // Read all the info
  VectorImageArray nu, ivkx, ivxk;
  ImageArray wkx, wxk;

  // Read all the warps nu (template to image i)
  for (int i = 0; i < N; i++)
    {
    if (i != k)
      {
      vimg_read(argv[iarg + i], nu[i]);
      vimg_read(argv[iarg + (N-1) + i], ivkx[i]);
      img_read(argv[iarg + 2 * (N-1) + i], wkx[i]);
      vimg_read(argv[iarg + 3 * (N-1) + i], ivxk[i]);
      img_read(argv[iarg + 4 * (N-1) + i], wxk[i]);
      }
    }

  // Read all the warps vkx (image x to image k)
  for (int i = 0; i < N; i++)
    vimg_read(argv[iarg++], vkx[i]);

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
