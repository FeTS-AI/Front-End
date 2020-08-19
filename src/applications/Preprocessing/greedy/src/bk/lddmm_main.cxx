#include <iostream>
#include <cstdio>
#include "lddmm_common.h"
#include "lddmm_data.h"

#include <itkImageFileReader.h>
#include <itkGaussianInterpolateImageFunction.h>
#include <itkResampleImageFilter.h>
#include <itkIdentityTransform.h>
#include <itkShrinkImageFilter.h>

using namespace std;

int usage()
{
  printf("lddmm: Paul's LDDMM implementation\n");
  printf("Usage: \n");
  printf("  lddmm DIM [options] fixed.nii moving.nii\n");
  printf("  lddmm DIM --test test_id test_params\n");
  printf("Options: \n");
  return -1;
}

template <class TFloat, uint VDim>
class Resampler
{
public:
  typedef itk::Image<TFloat, VDim> ImageType;

  Resampler(ImageType *image, uint factor)
    {
    m_Filter = FilterType::New();
    m_Filter->SetInput(image);

    typedef itk::IdentityTransform<TFloat, VDim> TranType;
    typename TranType::Pointer tran = TranType::New();
    m_Filter->SetTransform(tran);

    m_Func = FuncType::New();
    double sigma[VDim];
    for(uint i = 0; i < VDim; i++)
      sigma[i] = factor * image->GetSpacing()[i] * 1.2;
    m_Func->SetParameters(sigma, 4.0);

    m_Filter->SetInterpolator(m_Func);

    m_Shrink = ShrinkType::New();
    m_Shrink->SetInput(image);
    m_Shrink->SetShrinkFactors(factor);
    m_Shrink->Update();

    m_Filter->UseReferenceImageOn();
    m_Filter->SetReferenceImage(m_Shrink->GetOutput());
    m_Filter->Update();
    }

  ImageType *GetResult()
    { return m_Filter->GetOutput(); }
private:
  typedef itk::ShrinkImageFilter<ImageType, ImageType> ShrinkType;
  typedef itk::ResampleImageFilter<ImageType, ImageType, TFloat> FilterType;
  typedef itk::GaussianInterpolateImageFunction<ImageType, TFloat> FuncType;
  typename ShrinkType::Pointer m_Shrink;
  typename FilterType::Pointer m_Filter;
  typename FuncType::Pointer m_Func;
};

template <class TFloat, uint VDim>
int run_test(int argc, char *argv[])
{
  typedef LDDMMData<TFloat, VDim> LDDMM;
  typedef typename LDDMM::ImageType ImageType;
  typedef typename LDDMM::VectorImageType VectorImageType;


  int test_id = atoi(argv[1]);
  if(test_id == 0)
    {
    printf(
      "TEST LISTING:\n"
      "  1: FFT Convolution test for vector fields\n"
      "      --test 1 invec.nii compare.nii\n"
      "  2: Test image warping given velocity field\n"
      "      --test 2 nt warp%%02d.nii input.nii compare.nii\n"
      "  3: Test gradient w.r.t. variation\n"
      "      --test 3 fixed.nii moving.nii nt variation.nii\n"
    ); 

    return 0;
    }

  if(test_id == 1)
    {
    if(argc < 4) 
      throw itk::ExceptionObject("Wrong number of arguments for test 1");

    // Read images 
    typename VectorImageType::Pointer src, comp;
    LDDMM::vimg_read(argv[2], src);
    LDDMM::vimg_read(argv[3], comp);

    // Initialize LDDMM
    typename ImageType::Pointer dummy_fix, dummy_mov;
    LDDMM::alloc_img(dummy_fix, src);
    LDDMM::alloc_img(dummy_mov, src);
    LDDMM lddmm;
    LDDMM::init(lddmm, dummy_fix, dummy_mov, 4, 0.01, 1.0, 1.0);
    LDDMM::img_write(lddmm.f_kernel_sq, "test01_kernel.nii");

    // Perform operation
    typedef LDDMMFFTInterface<TFloat, VDim> FFT;
    FFT fft(dummy_fix);
    fft.convolution_fft(src, lddmm.f_kernel_sq, false, src);

    // Save (temp)
    LDDMM::vimg_write(src, "test01_result.nii");

    // Compare to target
    LDDMM::vimg_subtract_in_place(src, comp);
    TFloat diff = LDDMM::vimg_euclidean_norm_sq(src);

    printf("Difference with expected result: %f\n", diff);
    return (diff > 0);
    }
  else if(test_id == 2)
    {
    if(argc < 6) 
      throw itk::ExceptionObject("Wrong number of arguments for test 1");

    // Get the number of time steps
    uint nt = atoi(argv[2]);

    // Read input and output images
    typename ImageType::Pointer src, comp;
    LDDMM::img_read(argv[4], src);
    LDDMM::img_read(argv[5], comp);

    // Create a problem
    LDDMM lddmm;
    LDDMM::init(lddmm, src, src, nt, 0.01, 1.0, 1.0);

    // Read velocity field
    LDDMM::vfield_read(nt, argv[3], lddmm.v);

    // Integrate the forward transform
    lddmm.compute_semi_lagrangean_a();
    lddmm.integrate_phi_t0();

    // Warp the image using phi_10
    typename ImageType::Pointer res;
    LDDMM::alloc_img(res, src);
    LDDMM::interp_img(src, lddmm.f[nt-1], res);

    // Write the result
    LDDMM::img_write(res, "test02_result.nii");

    // Compare to target
    LDDMM::img_subtract_in_place(res, comp);
    TFloat diff = LDDMM::img_euclidean_norm_sq(res);

    printf("Difference with expected result: %f\n", diff);
    return (diff > 0);
    }
  else if(test_id == 3)
    {
    if(argc < 6) 
      throw itk::ExceptionObject("Wrong number of arguments for test 1");

    // Read input and output images
    typename ImageType::Pointer fix, mov;
    LDDMM::img_read(argv[2], fix);
    LDDMM::img_read(argv[3], mov);

    // Number of time steps
    uint nt = atoi(argv[4]);

    // Create a problem
    LDDMM p;
    LDDMM::init(p, fix, mov, nt, 0.01, 1.0, 1.0);

    // Iterate for 4 iterations (so that the starting point isn't zero
    LDDMMImageMatchingObjective<TFloat, VDim> obj(p);
    double ts = 0.01;
    for(uint k = 0; k < 4; k++)
      {
      TFloat fx = obj.compute_objective_and_gradient(p);
      for(uint it = 0; it < p.nt; it++)
        LDDMM::vimg_add_scaled_in_place(p.v[it], p.a[it], -ts);
      printf("Iter %04d     Obj: %10.10f\n", k, fx);
      }

    // Call this one more time, so we have v in p.v and dv in p.a
    TFloat fx = obj.compute_objective_and_gradient(p);

    // Read the random variation
    typename LDDMM::VelocityField h;
    LDDMM::alloc_vf(h, p.nt, p.fix);
    LDDMM::vfield_read(nt, argv[5], h);

    // gateaux_analytic = lddmm_vector_field_dot_product(dedvx, dedvy, varx, vary, p);
    
    // for expeps = -6:2
        
    //    eps = 10^expeps;
    
    //    E1 = lddmm_objective_and_gradient(vx - eps * varx, vy - eps * vary, p);
    //    E2 = lddmm_objective_and_gradient(vx + eps * varx, vy + eps * vary, p);

    //    gateaux_numeric = (E2 - E1) / (2 * eps);

    //    fprintf('Iter: %4i    Eps: %8d    Analytic: %12d     Numeric: %12d\n', ...
    //        i, eps, gateaux_analytic, gateaux_numeric);
    }
  else
    throw itk::ExceptionObject("Unknown test ID");
}

template <class TFloat, uint VDim>
int my_main(int argc, char *argv[])
{
  // Parse options, look for test option
  for(int i = 1; i < argc-1; i++)
    {
    if(!strcmp(argv[i], "--test"))
      return run_test<TFloat,VDim>(argc-i, argv+i);
    }

  // Normal processing
  char *fnmov = argv[argc-1];
  char *fnfix = argv[argc-2]; 
  // const char *outname = "lddmm.nii.gz";

  // Number of iterations
  uint n_res = 1;
  uint n_iter[] = {100};

  // Time steps
  uint nt = 10;

  // Read the images
  typedef LDDMMData<TFloat, VDim> LDDMM;
  typedef typename LDDMM::ImageType ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer readfix = ReaderType::New();
  readfix->SetFileName(fnfix);
  readfix->Update();
  
  typename ReaderType::Pointer readmov = ReaderType::New();
  readmov->SetFileName(fnmov);
  readmov->Update();
  
  // Loop over the multi-resolution levels
  for(int ires = n_res-1; ires >=0; ires--)
    {
    // Interpolate the input image and output image
    typename ImageType::Pointer ifix, imov;
    if(ires > 0) 
      {   
      Resampler<TFloat,VDim> rfix(readfix->GetOutput(), 1 << ires);
      Resampler<TFloat,VDim> rmov(readmov->GetOutput(), 1 << ires);
      ifix = rfix.GetResult();
      imov = rmov.GetResult();
      }
    else
      {
      ifix = readfix->GetOutput();
      imov = readmov->GetOutput();
      }

    // Create the problem
    LDDMM p;

    double avgdim = pow(ifix->GetBufferedRegion().GetNumberOfPixels() * 1.0, 1.0 / VDim);

    // LDDMM::init(p, ifix, imov, nt, 0.01, 1, 0.008);
    LDDMM::init(p, ifix, imov, nt, 0.01, 1, 1.0 / avgdim);

    // TODO: initialize with earlier set of velocity fields
    
    // Create an optimization problem
    LDDMMImageMatchingObjective<TFloat, VDim> obj(p);

    // Set time step
    double ts = 0.01;

    // Iterate
    for(uint k = 0; k < n_iter[ires]; k++)
      {
      // Compute objective and gradient. This will put the gradient of the function
      // into the 'a' array, retain the velocity field in 'v', and retain the inverse
      // transform (ft0) in the 'f' array. 
      TFloat fx = obj.compute_objective_and_gradient(p);

      // Update the solution by the time step (v = v - t * dv)
      for(uint it = 0; it < p.nt; it++)
        LDDMM::vimg_add_scaled_in_place(p.v[it], p.a[it], -ts);

      // Save the current image
      char fn[1024]; sprintf(fn, "lddmm_result_%04d.nii.gz", k);
      LDDMM::img_write(obj.Jt0, fn);
      
      // Print out the current status
      printf("Iter %04d     Obj: %10.10f\n", k, fx);
      }
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
