#include <iostream>
#include <cstdio>
#include "lddmm_common.h"
#include "lddmm_data.h"
#include "itkMultiThreader.h"

int usage()
{
  printf("macf_objective_and_delta: Paul's MACF implementation\n");
  printf("Usage: \n");
  printf("  macf_objgrad DIM N phiinv_i [phiinv_j] [psi_ij] [w_ij] [delta_i]\n");
  printf("Parameters: \n");
  printf("  phiinv_i :  A warp (from subject to template)\n");
  printf("  phiinv_j :  A list of N-1 warps (from other subjects to template)\n");
  printf("  psi_ij   :  A list of N-1 warps (from subject to subject)\n");
  printf("  w_ij     :  A list of N-1 weight maps\n");
  printf("  delta_i  :  Output image \n");
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

  printf("MultiThreader Info: MAX = %d   DEF = %d\n",
    itk::MultiThreader::GetGlobalMaximumNumberOfThreads(),
    itk::MultiThreader::GetGlobalDefaultNumberOfThreads());

  // No parameter checking!
  if(argc < 3) return usage();

  // Basic info
  unsigned int N = (unsigned int) atoi(argv[2]);
  if(argc < (int) (3 * N + 2)) return usage();

  // Get the lists of images
  char **fnphi = argv + 3;
  char **fnpsi = fnphi + N;
  char **fnw = fnpsi + N - 1;
  char **fndelta = fnw + N - 1;

  // Read i-th warp, store in delta
  VectorImagePointer delta;
  LDDMM::vimg_read(fnphi[0], delta);

  // Temporary velocity field for computations
  VectorImagePointer q;
  LDDMM::alloc_vimg(q, delta);

  // Iterate over all other warps
  for(unsigned int i = 0; i < N-1; i++)
    {
    VectorImagePointer psi_ij, phi_j;
    ImagePointer w_ij;

    // Progress report
    printf("Updating Delta_i by %s, %s, %s\n", fnphi[i+1], fnpsi[i], fnw[i]);

    // Read the data
    LDDMM::vimg_read(fnphi[i+1], phi_j);
    LDDMM::vimg_read(fnpsi[i], psi_ij);
    LDDMM::img_read(fnw[i], w_ij);

    // Compose transformations
    LDDMM::interp_vimg(phi_j, psi_ij, 1.0, true, q);
    LDDMM::vimg_add_in_place(q, psi_ij);

    // Scale by W
    LDDMM::vimg_multiply_in_place(q, w_ij);

    // Subtract from current delta estimate
    LDDMM::vimg_subtract_in_place(delta, q);
    }

  // Compute the energy for this iteration
  double energy = LDDMM::vimg_euclidean_norm_sq(delta);
  printf("ENERGY = %lf\n", energy);

  // Write the delta
  LDDMM::vimg_write(delta, fndelta[0]);

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
    my_main<float, 2>(argc, argv);
  else if(ndim == 3)
    my_main<float, 3>(argc, argv);
  else
    return usage();

}
