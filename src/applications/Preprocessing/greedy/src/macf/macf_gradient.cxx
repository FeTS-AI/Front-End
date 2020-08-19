#include <iostream>
#include <cstdio>
#include "lddmm_common.h"
#include "lddmm_data.h"

int usage()
{
  printf("macf_gradient: Paul's MACF implementation\n");
  printf("Usage: \n");
  printf("  macf_gradient DIM N sigma dt delta_i [delta_j] [psi_ji] [J_psi_ji] [w_ji] phi_i phi_i_out\n");
  printf("Parameters: \n"); 
  printf("  DIM      :  Image dimensionality\n");        // 1
  printf("  N        :  Number of images\n");        // 1
  printf("  sigma    :  Gaussian smoothing (VOXEL units)\n");        // 1
  printf("  dt       :  Step size\n"); // 1
  printf("  delta_i  :  A delta warp (from macf_objective_and_delta )\n"); // 1
  printf("  delta_j  :  A list of N-1 deltas (for other subjects)\n"); // N-1
  printf("  psi_ji   :  A list of N-1 warps (from subject to subject)\n");
  printf("  J_psi_ji :  A list of N-1 jacobian determinants of above warps\n");
  printf("  w_ji     :  A list of N-1 weight maps\n");
  printf("  phi_i    :  Mapping from subject to template from last step\n"); // 1
  printf("  phi_out_i:  Output image\n"); // 1
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
  if(argc < 3) return usage();

  // Basic info
  unsigned int N = (unsigned int) atoi(argv[2]);
  if(argc < (int) (4 * N + 4)) return usage();

  // Simple parameters
  double sigma = atof(argv[3]);
  double dt = atof(argv[4]);

  // Get the lists of images
  char **fndelta = argv + 5;
  char **fnpsi = fndelta + N;
  char **fnJpsi = fnpsi + N - 1;
  char **fnw = fnJpsi + N - 1;
  char **fnphi = fnw + N - 1;

  // Read i-th warp, store in delta
  VectorImagePointer grad;
  LDDMM::vimg_read(fndelta[0], grad);

  // Temporary velocity field for computations
  VectorImagePointer q;
  LDDMM::alloc_vimg(q, grad);

  // Iterate over all other warps
  for(unsigned int j = 0; j < N-1; j++)
    {
    VectorImagePointer psi_ji, delta_j;
    ImagePointer w_ji, Jpsi_ji;

    // Progress report
    printf("Updating Grad_i by %s, %s, %s\n", fndelta[j+1], fnpsi[j], fnw[j]);

    // Read the data
    LDDMM::vimg_read(fndelta[j+1], delta_j);
    LDDMM::vimg_read(fnpsi[j], psi_ji);
    LDDMM::img_read(fnw[j], w_ji);
    LDDMM::img_read(fnJpsi[j], Jpsi_ji);

    // Scale delta_j by w
    LDDMM::vimg_multiply_in_place(delta_j, w_ji);

    // Interpolate using psi_ji
    LDDMM::interp_vimg(delta_j, psi_ji, 1.0, true, q);

    // Scale by the Jacobian of this map
    LDDMM::vimg_multiply_in_place(q, Jpsi_ji);

    // Subtract from grad
    LDDMM::vimg_subtract_in_place(grad, q);
    }

  // Smooth the gradient and scale by timestep
  LDDMM::vimg_gaussian_smooth(grad, false, sigma, q);
  LDDMM::vimg_scale_in_place(q, -dt);

  // Read the current transformation
  VectorImagePointer phi_i;
  LDDMM::vimg_read(fnphi[0], phi_i); 

  // Compose with the last phi to get the new phi
  VectorImagePointer phiout = grad;
  LDDMM::interp_vimg(phi_i, q, 1.0, true, phiout);
  LDDMM::vimg_add_in_place(phiout, q);

  // Save the new phi
  LDDMM::vimg_write(phiout, fnphi[1]);

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
