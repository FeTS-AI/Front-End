#include <iostream>
#include <cstdio>
#include <vector>
#include <string>
#include <algorithm>

#include "lddmm_common.h"
#include "lddmm_data.h"

#include <itkImageFileReader.h>
#include <itkGaussianInterpolateImageFunction.h>
#include <itkResampleImageFilter.h>
#include <itkIdentityTransform.h>
#include <itkShrinkImageFilter.h>



int usage()
{
  printf("greedy: Paul's greedy diffeomorphic registration implementation\n");
  printf("Usage: \n");
  printf("  greedy [options]\n");
  printf("Required options: \n");
  printf("  -d DIM                      : Number of image dimensions\n");
  printf("  -i fixed.nii moving.nii     : Image pair (may be repeated)\n");
  printf("  -o output.nii               : Output warp file\n");
  printf("Optional: \n");
  printf("  -w weight                   : weight of the next -i pair\n");
  printf("  -e epsilon                  : step size (default = 1.0)\n");
  printf("  -s sigma                    : smoothing for the greedy update step (3.0)\n");
  printf("  -n number                   : number of iterations (200) \n");
  printf("  -dump-moving                : dump moving image at each iter\n");
  printf("  -dump-freq N                : dump frequency\n");
  return -1;
}

struct ImagePairSpec
{
  std::string fixed;
  std::string moving;
  double weight;
};

struct GreedyParameters
{
  std::vector<ImagePairSpec> inputs;
  std::string output;
  unsigned int dim; 

  bool flag_dump_moving;
  int dump_frequency;
  double epsilon, sigma;
  int niter;
};


template <unsigned int VDim, typename TReal = double>
class GreedyApproach
{
public:

  typedef LDDMMData<TReal, VDim> LDDMMType;
  typedef typename LDDMMType::ImageType ImageType;
  typedef typename LDDMMType::ImagePointer ImagePointer;
  typedef typename LDDMMType::VectorImageType VectorImageType;
  typedef typename LDDMMType::VectorImagePointer VectorImagePointer;

  struct ImagePair {
    ImagePointer fixed, moving;
    VectorImagePointer grad_moving;
    double weight;
  };

  static int Run(GreedyParameters &param);

protected:
};

/**
 * This is the main function of the GreedyApproach algorithm
 */
template <unsigned int VDim, typename TReal>
int GreedyApproach<VDim, TReal>
::Run(GreedyParameters &param)
{
  // Image pairs to register
  std::vector<ImagePair> img;

  // Read the input images
  for(int i = 0; i < param.inputs.size(); i++)
    {
    ImagePair ip;

    // Read fixed
    typedef itk::ImageFileReader<ImageType> ReaderType;
    typename ReaderType::Pointer readfix = ReaderType::New();
    readfix->SetFileName(param.inputs[i].fixed);
    readfix->Update();
    ip.fixed = readfix->GetOutput();

    // Read moving
    typedef itk::ImageFileReader<ImageType> ReaderType;
    typename ReaderType::Pointer readmov = ReaderType::New();
    readmov->SetFileName(param.inputs[i].moving);
    readmov->Update();
    ip.moving = readmov->GetOutput();

    // Allocate the gradient
    LDDMMType::alloc_vimg(ip.grad_moving, ip.moving);

    // Precompute the gradient of the moving images. There should be some
    // smoothing of the input images before applying this computation!
    LDDMMType::image_gradient(ip.moving, ip.grad_moving);

    // Set weight
    ip.weight = param.inputs[i].weight;

    // Append
    img.push_back(ip);
    }

  // Reference space
  ImagePointer refspace = img.front().fixed;

  // Initialize the displacement to identity
  VectorImagePointer uk = VectorImageType::New();
  LDDMMType::alloc_vimg(uk, refspace);

  // Initialize the intermediate data
  ImagePointer iTemp = ImageType::New();
  LDDMMType::alloc_img(iTemp, refspace);

  VectorImagePointer viTemp = VectorImageType::New();
  LDDMMType::alloc_vimg(viTemp, refspace);

  VectorImagePointer uk1 = VectorImageType::New();
  LDDMMType::alloc_vimg(uk1, refspace);

  // Iterate
  for(unsigned int iter = 0; iter < param.niter; iter++)
    {
    // Initialize u(k+1) to zero
    uk1->FillBuffer(typename LDDMMType::Vec(0.0));

    // Initialize the energy computation
    double total_energy = 0.0;

    // Add all the derivative terms
    for(int j = 0; j < img.size(); j++)
      {
      // Interpolate each moving image
      LDDMMType::interp_img(img[j].moving, uk, iTemp);

      // Dump the moving image?
      if(param.flag_dump_moving && 0 == iter % param.dump_frequency)
        {
        char fname[256];
        sprintf(fname, "dump_moving_%02d_iter%04d.nii.gz", j, iter);
        LDDMMType::img_write(iTemp, fname);
        }

      // Subtract the fixed image
      LDDMMType::img_subtract_in_place(iTemp, img[j].fixed);

      // Record the norm of the difference image
      total_energy += img[j].weight * LDDMMType::img_euclidean_norm_sq(iTemp);

      // Interpolate the gradient of the moving image
      LDDMMType::interp_vimg(img[j].grad_moving, uk, 1.0, viTemp);
      LDDMMType::vimg_multiply_in_place(viTemp, iTemp);

      // Accumulate to the force
      LDDMMType::vimg_add_scaled_in_place(uk1, viTemp, -img[j].weight * param.epsilon);
      }

    // We have now computed the gradient vector field. Next, we smooth it 
    LDDMMType::vimg_smooth(uk1, viTemp, param.sigma);

    // Write Uk1
    if(param.flag_dump_moving && 0 == iter % param.dump_frequency)
      {
      char fname[256];
      sprintf(fname, "dump_graduent_iter%04d.nii.gz", iter);
      LDDMMType::vimg_write(uk1, fname);
      sprintf(fname, "dump_optflow_iter%04d.nii.gz", iter);
      LDDMMType::vimg_write(viTemp, fname);
      }

    // Compute the updated deformation field - in uk1
    LDDMMType::interp_vimg(uk, viTemp, 1.0, uk1);
    LDDMMType::vimg_add_in_place(uk1, viTemp);

    if(param.flag_dump_moving && 0 == iter % param.dump_frequency)
      {
      char fname[256];
      sprintf(fname, "dump_uk1_iter%04d.nii.gz", iter);
      LDDMMType::vimg_write(uk1, fname);
      }

    // Swap uk and uk1 pointers
    VectorImagePointer tmpptr = uk1; uk1 = uk; uk = tmpptr;
    // LDDMMType::vimg_smooth(uk1, uk, 1.0);

    // Compute the Jacobian determinant of the updated field (temporary)
    LDDMMType::field_jacobian_det(uk, iTemp);
    TReal jac_min, jac_max;
    LDDMMType::img_min_max(iTemp, jac_min, jac_max);

    // Report the energy
    printf("Iter %5d:    Energy = %8.4f     DetJac Range: %8.4f  to %8.4f \n", iter, total_energy, jac_min, jac_max);
    }

  // Write the resulting transformation field
  LDDMMType::vimg_write(uk, param.output.c_str());

  return 0;
}



int main(int argc, char *argv[])
{
  GreedyParameters param;
  double current_weight = 1.0;

  param.dim = 2;
  param.flag_dump_moving = false;
  param.dump_frequency = 1;
  param.epsilon = 1.0;
  param.sigma = 3.0;
  param.niter = 200;

  if(argc < 3)
    return usage();

  for(int i = 1; i < argc; ++i)
    {
    std::string arg = argv[i];
    if(arg == "-d")
      {
      param.dim = atoi(argv[++i]);
      }
    else if(arg == "-n")
      {
      param.niter = atoi(argv[++i]);
      }
    else if(arg == "-w")
      {
      current_weight = atof(argv[++i]);
      }
    else if(arg == "-e")
      {
      param.epsilon = atof(argv[++i]);
      }
    else if(arg == "-s")
      {
      param.sigma = atof(argv[++i]);
      }
    else if(arg == "-i")
      {
      ImagePairSpec ip;
      ip.weight = current_weight;
      ip.fixed = argv[++i];
      ip.moving = argv[++i];
      param.inputs.push_back(ip);
      }
    else if(arg == "-o")
      {
      param.output = argv[++i];
      }
    else if(arg == "-dump-moving")
      {
      param.flag_dump_moving = true;
      }
    else if(arg == "-dump-frequency")
      {
      param.dump_frequency = atoi(argv[++i]);
      }
    else
      {
      std::cerr << "Unknown parameter " << arg << std::endl;
      return -1;
      }
    }

  switch(param.dim)
    {
    case 2: return GreedyApproach<2>::Run(param); break;
    case 3: return GreedyApproach<3>::Run(param); break;
    case 4: return GreedyApproach<4>::Run(param); break;
    default: 
          std::cerr << "Wrong dimensionality" << std::endl;
          return -1;
    }
}
