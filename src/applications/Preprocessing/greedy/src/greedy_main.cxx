/*=========================================================================

  Program:   ALFABIS fast medical image registration programs
  Language:  C++
  Website:   github.com/pyushkevich/greedy
  Copyright (c) Paul Yushkevich, University of Pennsylvania. All rights reserved.

  This program is part of ALFABIS: Adaptive Large-Scale Framework for
  Automatic Biomedical Image Segmentation.

  ALFABIS development is funded by the NIH grant R01 EB017255.

  ALFABIS is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ALFABIS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ALFABIS.  If not, see <http://www.gnu.org/licenses/>.

=========================================================================*/

#include "GreedyAPI.h"
#include "CommandLineHelper.h"
#include "cbicaCmdParser.h"

#include <iostream>
#include <sstream>
#include <cstdio>
#include <vector>
#include <string>
#include <algorithm>
#include <cerrno>

#include "lddmm_common.h"
#include "lddmm_data.h"

#include <itkImageFileReader.h>
#include <itkAffineTransform.h>
#include <itkTransformFactory.h>
#include <itkTimeProbe.h>

#include "MultiImageRegistrationHelper.h"
#include "FastWarpCompositeImageFilter.h"
#include <vnl/vnl_cost_function.h>
#include <vnl/vnl_random.h>
#include <vnl/algo/vnl_powell.h>
#include <vnl/algo/vnl_svd.h>
#include <vnl/vnl_trace.h>

extern const char *GreedyVersionInfo;

int usage()
{
  printf("greedy: Paul's greedy diffeomorphic registration implementation\n");
  printf("Usage: \n");
  printf("  greedy [options]\n");
  printf("Required options: \n");
  printf("  -d DIM                 : Number of image dimensions\n");
  printf("  -i fix.nii mov.nii     : Image pair (may be repeated)\n");
  printf("  -o output.nii          : Output file\n");
  printf("Mode specification: \n");
  printf("  -a                     : Perform affine registration and save to output (-o)\n");
  printf("  -brute radius          : Perform a brute force search around each voxel \n");
  printf("  -moments <1|2>         : Perform moments of inertia rigid alignment of given order.\n");
  printf("                               order 1 matches center of mass only\n");
  printf("                               order 2 matches second-order moments of inertia tensors\n");
  printf("  -r [tran_spec]         : Reslice images instead of doing registration \n");
  printf("                               tran_spec is a series of warps, affine matrices\n");
  printf("  -iw inwarp outwarp     : Invert previously computed warp\n");
  printf("  -root inwarp outwarp N : Convert 2^N-th root of a warp \n");
  printf("  -jac inwarp outjac     : Compute the Jacobian determinant of the warp \n");
  printf("Options in deformable / affine mode: \n");
  printf("  -w weight              : weight of the next -i pair\n");
  printf("  -m metric              : metric for the entire registration\n");
  printf("                               SSD:          sum of square differences (default)\n");
  printf("                               MI:           mutual information\n");
  printf("                               NMI:          normalized mutual information\n");
  printf("                               NCC <radius>: normalized cross-correlation\n");
  printf("                               MAHAL:        Mahalanobis distance to target warp\n");
  printf("  -e epsilon             : step size (default = 1.0), \n");
  printf("                               may also be specified per level (e.g. 0.3x0.1)\n");
  printf("  -n NxNxN               : number of iterations per level of multi-res (100x100) \n");
  printf("  -threads N             : set the number of allowed concurrent threads\n");
  printf("  -gm mask.nii           : mask for gradient computation\n");
  printf("  -gm-trim <radius>      : generate mask for gradient computation by trimming the extent\n");
  printf("                           of the fixed image by given radius. This is useful during affine\n");
  printf("                           registration with the NCC metric when the background of your images\n");
  printf("                           is non-zero. The radius should match that of the NCC metric.");
  printf("  -mm mask.nii           : mask for the moving image\n");
  printf("  -it filenames          : sequence of transforms to apply to the moving image first \n");
  printf("Specific to deformable mode: \n");
  printf("  -tscale MODE           : time step behavior mode: CONST, SCALE [def], SCALEDOWN\n");
  printf("  -s sigma1 sigma2       : smoothing for the greedy update step. Must specify units,\n");
  printf("                           either `vox` or `mm`. Default: 1.732vox, 0.7071vox\n");
  printf("  -oinv image.nii        : compute and write the inverse of the warp field into image.nii\n");
  printf("  -oroot image.nii       : compute and write the (2^N-th) root of the warp field into image.nii, where\n");
  printf("                           N is the value of the -exp option. In stational velocity mode, it is advised\n");
  printf("                           to output the root warp, since it is used internally to represent the deformation\n");
  printf("  -wp VALUE              : Saved warp precision (in voxels; def=0.1; 0 for no compression).\n");
  printf("  -noise VALUE           : Standard deviation of white noise added to moving/fixed images when \n");
  printf("                           using NCC metric. Relative to intensity range. Def=0.001\n");
  printf("  -exp N                 : The exponent used for warp inversion, root computation, and in stationary \n");
  printf("                           velocity field (Diff Demons) mode. N is a positive integer (default = 6) \n");
  printf("  -sv                    : Performs registration using the stationary velocity model, similar to diffeomoprhic \n");
  printf("                           Demons (Vercauteren 2008 MICCAI). Internally, the deformation field is \n");
  printf("                           represented as 2^N self-compositions of a small deformation and \n");
  printf("                           greedy updates are applied to this deformation. N is specified with the -exp \n");
  printf("                           option (6 is a good number). This mode results in better behaved\n");
  printf("                           deformation fields and Jacobians than the pure greedy approach.\n");
  printf("  -svlb                  : Same as -sv but uses the more accurate but also more expensive \n");
  printf("                           update of v, v <- v + u + [v,u]. Experimental feature \n");
  printf("  -id image.nii          : Specifies the initial warp to start iteration from. In stationary mode, this \n");
  printf("                           is the initial stationary velocity field (output by -oroot option)\n");
  printf("Initial transform specification: \n");
  printf("  -ia filename           : initial affine matrix for optimization (not the same as -it) \n");
  printf("  -ia-identity           : initialize affine matrix based on NIFTI headers \n");
  printf("  -ia-image-centers      : initialize affine matrix based on matching image centers \n");
  printf("  -ia-image-side CODE    : initialize affine matrix based on matching center of one image side \n");
  printf("  -ia-moments <1|2>      : initialize affine matrix based on matching moments of inertia\n");
  printf("Specific to affine mode (-a):\n");
  printf("  -dof N                 : Degrees of freedom for affine reg. 6=rigid, 12=affine\n");
  printf("  -jitter sigma          : Jitter (in voxel units) applied to sample points (def: 0.5)\n");
  printf("  -search N <rot> <tran> : Random search over rigid transforms (N iter) before starting optimization\n");
  printf("                           'rot' may be the standard deviation of the random rotation angle (degrees) or \n");
  printf("                           keyword 'any' (any rotation) or 'flip' (any rotation or flip). \n");
  printf("                           'tran' is the standard deviation of the random offset, in physical units. \n");
  printf("Specific to moments of inertia mode (-moments 2): \n");
  printf("  -det <-1|1>            : Force the determinant of transform to be either 1 (no flip) or -1 (flip)\n");
  printf("  -cov-id                : Assume identity covariance (match centers and do flips only, no rotation)\n");
  printf("Specific to reslice mode (-r): \n");
  printf("  -rf fixed.nii          : fixed image for reslicing\n");
  printf("  -rm mov.nii out.nii    : moving/output image pair (may be repeated)\n");
  printf("  -rs mov.vtk out.vtk    : moving/output surface pair (vertices are warped from fixed space to moving)\n");
  printf("  -ri interp_mode        : interpolation for the next pair (NN, LINEAR*, LABEL sigma)\n");
  printf("  -rb value              : background (i.e. outside) intensity for the next pair (default 0)\n");
  printf("  -rc outwarp            : write composed transforms to outwarp \n");
  printf("  -rj outjacobian        : write Jacobian determinant image to outjacobian \n");
  printf("For developers: \n");
  printf("  -debug-deriv           : enable periodic checks of derivatives (debug) \n");
  printf("  -debug-deriv-eps       : epsilon for derivative debugging \n");
  printf("  -debug-aff-obj         : plot affine objective in neighborhood of -ia matrix \n");
  printf("  -dump-moving           : dump moving image at each iter\n");
  printf("  -dump-freq N           : dump frequency\n");
  printf("  -powell                : use Powell's method instead of LGBFS\n");
  printf("  -float                 : use single precision floating point (off by default)\n");
  printf("  -version               : print version info\n");

  return -1;
}



template <unsigned int VDim, typename TReal>
class GreedyRunner
{
public:
  static int Run(GreedyParameters &param)
  {
    // Use the threads parameter
    GreedyApproach<VDim, TReal> greedy;   
    return greedy.Run(param);
  }
};





int main(int argc, char *argv[])
{
  GreedyParameters param;
  GreedyParameters::SetToDefaults(param);

  if(argc < 2)
    return usage();

  try
  {
    CommandLineHelper cl(argc, argv);
    while(!cl.is_at_end())
      {
      // Read the next command
      std::string cmd = cl.read_command();

      // Parse generic commands
      if(cmd == "-version")
        {
        std::cout << GreedyVersionInfo << std::endl;
        exit(0);
        }
      else if(cmd == "-h" || cmd == "-help" || cmd == "--help")
        {
        return usage();
        }
      else if(!param.ParseCommandLine(cmd, cl))
        {
        std::cerr << "Unknown parameter " << cmd << std::endl;
        return -1;
        }
      }

    // Some parameters may be specified as either vector or scalar, and need to be verified
    if(!param.epsilon_per_level.CheckSize(param.iter_per_level.size()))
       throw GreedyException("Mismatch in size of vectors supplied with -n and -e options");

    // Run the main code
    if(param.flag_float_math)
      {
      switch(param.dim)
        {
        case 2: return GreedyRunner<2, float>::Run(param); break;
        case 3: return GreedyRunner<3, float>::Run(param); break;
        case 4: return GreedyRunner<4, float>::Run(param); break;
        default: throw GreedyException("Wrong number of dimensions requested: %d", param.dim);
        }
      }
    else
      {
      switch(param.dim)
        {
        case 2: return GreedyRunner<2, double>::Run(param); break;
        case 3: return GreedyRunner<3, double>::Run(param); break;
        case 4: return GreedyRunner<4, double>::Run(param); break;
        default: throw GreedyException("Wrong number of dimensions requested: %d", param.dim);
        }
      }
  }
  catch(std::exception &exc)
  {
    std::cerr << "ABORTING PROGRAM DUE TO RUNTIME EXCEPTION -- "
              << exc.what() << std::endl;
    return -1;
  }
}
