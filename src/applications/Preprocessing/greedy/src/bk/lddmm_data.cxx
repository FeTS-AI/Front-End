#include "lddmm_data.h"
#include "itkImageRegionIterator.h"
#include "SimpleWarpImageFilter.h"
#include "itkNumericTraitsCovariantVectorPixel.h"
#include "itkOptVectorLinearInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkAddImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkGradientImageFilter.h"
#include "itkUnaryFunctorImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkVectorImage.h"
#include "itkDisplacementFieldJacobianDeterminantFilter.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkMinimumMaximumImageFilter.h"

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::alloc_vf(VelocityField &vf, uint nt, ImageBaseType *ref)
{
  vf.resize(nt);
  for(uint i = 0; i < nt; i++)
    alloc_vimg(vf[i], ref);
}

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::alloc_vimg(VectorImagePointer &img, ImageBaseType *ref)
{
  img = VectorImageType::New();
  img->SetRegions(ref->GetBufferedRegion());
  img->CopyInformation(ref);
  img->Allocate();
  img->FillBuffer(Vec(0.0));
}

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::alloc_img(ImagePointer &img, ImageBaseType *ref)
{
  img = ImageType::New();
  img->SetRegions(ref->GetBufferedRegion());
  img->CopyInformation(ref);
  img->Allocate();
  img->FillBuffer(0.0);
}

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::init(LDDMMData<TFloat, VDim> &p, 
  ImageType *fix, ImageType *mov, 
  uint nt, double alpha, double gamma, double sigma)
{
  p.fix = fix;
  p.mov = mov;
  p.alpha = alpha;
  p.sigma = sigma;
  p.gamma = gamma;
  p.nt = nt;
  p.dt = 1.0 / (nt - 1.0);
  p.sigma_sq = sigma * sigma;

  // Initialize N and R
  p.r = fix->GetBufferedRegion();
  p.nv = fix->GetBufferedRegion().GetNumberOfPixels();
  for(uint i = 0; i < VDim; i++)
    p.n[i] = p.r.GetSize()[i];

  // Initialize the velocity fields
  alloc_vf(p.v, nt, fix);
  alloc_vf(p.a, nt, fix);
  alloc_vf(p.f, nt, fix);

  // Initialize kernel terms
  alloc_img(p.f_kernel, fix);
  alloc_img(p.f_kernel_sq, fix);

  // Compute these images
  ImageIterator it(p.f_kernel, p.r), itsq(p.f_kernel_sq, p.r);
  for(; !it.IsAtEnd(); ++it, ++itsq)
    {
    TFloat val = 0.0;
    for(uint j = 0; j < VDim; j++)
      val += 1.0 - cos(it.GetIndex()[j] * 2.0 * vnl_math::pi / p.n[j]);
    it.Set(p.gamma + 2.0 * p.alpha * p.nv * val);
    itsq.Set(it.Get() * it.Get());
    }

  // Initialize temporary vector field
  alloc_vimg(p.vtmp, fix);
}

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::interp_vimg(
  VectorImageType *data, VectorImageType *field, 
  TFloat def_scale, VectorImageType *out)
{
  // Create a warp filter
  typedef itk::SimpleWarpImageFilter<
    VectorImageType, VectorImageType, VectorImageType, TFloat> WarpFilterType;
  typename WarpFilterType::Pointer flt = WarpFilterType::New();

  // Create an interpolation function
  typedef itk::OptVectorLinearInterpolateImageFunction<
    VectorImageType, TFloat> InterpType;
  typename InterpType::Pointer func = InterpType::New();

  // Graft output of the warp filter
  flt->GraftOutput(out);

  // Set inputs
  flt->SetInput(data);
  flt->SetInterpolator(func.GetPointer());
  flt->SetDeformationField(field);
  flt->SetDeformationScaling(def_scale);
  flt->Update();

}

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::interp_img(ImageType *data, VectorImageType *field, ImageType *out)
{
  // Create a warp filter
  typedef itk::SimpleWarpImageFilter<
    ImageType, ImageType, VectorImageType, TFloat> WarpFilterType;
  typename WarpFilterType::Pointer flt = WarpFilterType::New();

  // Create an interpolation function
  typedef itk::LinearInterpolateImageFunction<ImageType, TFloat> InterpType;
  typename InterpType::Pointer func = InterpType::New();

  // Graft output of the warp filter
  flt->GraftOutput(out);

  // Set inputs
  flt->SetInput(data);
  flt->SetInterpolator(func.GetPointer());
  flt->SetDeformationField(field);
  flt->SetDeformationScaling(1.0);
  flt->Update();
}

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::vimg_add_in_place(VectorImagePointer &trg, VectorImageType *a)
{
  typedef itk::AddImageFilter<VectorImageType> AddFilter;
  typename AddFilter::Pointer flt = AddFilter::New();
  flt->SetInput(0,trg);
  flt->SetInput(1,a);
  flt->InPlaceOn();
  flt->Update();
  trg = flt->GetOutput();
}

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::vimg_subtract_in_place(VectorImagePointer &trg, VectorImageType *a)
{
  typedef itk::SubtractImageFilter<VectorImageType> SubtractFilter;
  typename SubtractFilter::Pointer flt = SubtractFilter::New();
  flt->SetInput(0,trg);
  flt->SetInput(1,a);
  flt->InPlaceOn();
  flt->Update();
  trg = flt->GetOutput();
}

// Scalar math
template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::vimg_multiply_in_place(VectorImagePointer &trg, ImageType *s)
{
  typedef itk::MultiplyImageFilter<
    VectorImageType, ImageType, VectorImageType> MultiplyFilter;
  typename MultiplyFilter::Pointer flt = MultiplyFilter::New();
  flt->SetInput1(trg);
  flt->SetInput2(s);
  flt->InPlaceOn();
  flt->Update();
  trg = flt->GetOutput();
}

// Scalar math
template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::img_subtract_in_place(ImagePointer &trg, ImageType *a)
{
  typedef itk::SubtractImageFilter<ImageType> SubtractFilter;
  typename SubtractFilter::Pointer flt = SubtractFilter::New();
  flt->SetInput(0,trg);
  flt->SetInput(1,a);
  flt->InPlaceOn();
  flt->Update();
  trg = flt->GetOutput();
}

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::img_multiply_in_place(ImagePointer &trg, ImageType *a)
{
  typedef itk::MultiplyImageFilter<ImageType> MultiplyFilter;
  typename MultiplyFilter::Pointer flt = MultiplyFilter::New();
  flt->SetInput(0,trg);
  flt->SetInput(1,a);
  flt->InPlaceOn();
  flt->Update();
  trg = flt->GetOutput();
}

template <class TFloat, uint VDim>
TFloat 
LDDMMData<TFloat, VDim>
::vimg_euclidean_norm_sq(VectorImageType *trg)
{
  // Add all voxels in the image
  double accum = 0.0;
  typedef itk::ImageRegionIterator<VectorImageType> Iter;
  for(Iter it(trg, trg->GetBufferedRegion()); !it.IsAtEnd(); ++it)
    {
    for(uint d = 0; d < VDim; d++)
      accum += it.Value()[d] * it.Value()[d];
    }
  return (TFloat) accum;
}

template <class TFloat, uint VDim>
TFloat 
LDDMMData<TFloat, VDim>
::img_euclidean_norm_sq(ImageType *trg)
{
  // Add all voxels in the image
  double accum = 0.0;
  typedef itk::ImageRegionIterator<ImageType> Iter;
  for(Iter it(trg, trg->GetBufferedRegion()); !it.IsAtEnd(); ++it)
    { accum += it.Value() * it.Value(); }
  return (TFloat) accum;
}

template <class TFloat, uint VDim>
TFloat 
LDDMMData<TFloat, VDim>
::img_voxel_sum(ImageType *trg)
{
  // Add all voxels in the image
  double accum = 0.0;
  typedef itk::ImageRegionIterator<ImageType> Iter;
  for(Iter it(trg, trg->GetBufferedRegion()); !it.IsAtEnd(); ++it)
    accum += it.Value();
  return (TFloat) accum;
}

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::img_min_max(ImageType *src, TFloat &out_min, TFloat &out_max)
{
  // Add all voxels in the image
  typedef itk::MinimumMaximumImageFilter<ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput(src);
  filter->Update();
  out_min = filter->GetMinimum();
  out_max = filter->GetMaximum();
}


template <class TFloat, uint VDim>
class VectorScaleFunctor
{
public:
  VectorScaleFunctor() { this->Scale = 1.0; }
  typedef itk::CovariantVector<TFloat,VDim> Vec;

  Vec operator() (const Vec &x)
    { return x * Scale; }

  bool operator== (const VectorScaleFunctor<TFloat, VDim> &other)
    { return Scale == other.Scale; }

  bool operator!= (const VectorScaleFunctor<TFloat, VDim> &other)
    { return Scale != other.Scale; }

  TFloat Scale;
};

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::vimg_scale_in_place(VectorImagePointer &trg, TFloat s)
{
  typedef VectorScaleFunctor<TFloat, VDim> Functor;
  typedef itk::UnaryFunctorImageFilter<
    VectorImageType, VectorImageType, Functor> Filter;
  typename Filter::Pointer flt = Filter::New();

  Functor func;
  func.Scale = s;
  flt->SetFunctor(func);
  flt->SetInput(trg);
  flt->InPlaceOn();
  flt->Update();

  trg = flt->GetOutput();
}

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::vimg_scale(VectorImageType*src, TFloat s, VectorImagePointer &trg)
{
  typedef VectorScaleFunctor<TFloat, VDim> Functor;
  typedef itk::UnaryFunctorImageFilter<
    VectorImageType, VectorImageType, Functor> Filter;
  typename Filter::Pointer flt = Filter::New();

  Functor func;
  func.Scale = s;
  flt->SetFunctor(func);
  flt->SetInput(src);
  flt->InPlaceOff();
  flt->GraftOutput(trg);
  flt->Update();
}

template <class TFloat, uint VDim>
class VectorScaleAddFunctor
{
public:
  VectorScaleAddFunctor() { this->Scale = 1.0; }
  typedef itk::CovariantVector<TFloat,VDim> Vec;

  Vec operator() (const Vec &x, const Vec &y)
    { return x + y * Scale; }

  bool operator== (const VectorScaleAddFunctor<TFloat, VDim> &other)
    { return Scale == other.Scale; }

  bool operator!= (const VectorScaleAddFunctor<TFloat, VDim> &other)
    { return Scale != other.Scale; }

  TFloat Scale;
};

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::vimg_add_scaled_in_place(VectorImagePointer &trg, VectorImageType *a, TFloat s)
{
  typedef VectorScaleAddFunctor<TFloat, VDim> Functor;
  typedef itk::BinaryFunctorImageFilter<
    VectorImageType, VectorImageType, VectorImageType, Functor> Filter;
  typename Filter::Pointer flt = Filter::New();

  Functor func;
  func.Scale = s;
  flt->SetFunctor(func);
  flt->SetInput1(trg);
  flt->SetInput2(a);
  flt->InPlaceOn();
  flt->Update();

  trg = flt->GetOutput();
}

template <class TFloat, uint VDim>
class VectorDotProduct
{
public:
  VectorDotProduct() {}
  typedef itk::CovariantVector<TFloat,VDim> Vec;

  TFloat operator() (const Vec &x, const Vec &y)
    {
    TFloat dp = 0.0;
    for(uint d = 0; d < VDim; d++)
      dp += x[d] * y[d];
    return dp;
    }

  bool operator== (const VectorDotProduct<TFloat, VDim> &other)
    { return true; }

  bool operator!= (const VectorDotProduct<TFloat, VDim> &other)
    { return false; }
};

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::vimg_euclidean_inner_product(ImagePointer &trg, VectorImageType *a, VectorImageType *b)
{
  typedef VectorDotProduct<TFloat, VDim> Functor;
  typedef itk::BinaryFunctorImageFilter<
    VectorImageType, VectorImageType, ImageType, Functor> Filter;
  typename Filter::Pointer flt = Filter::New();

  Functor func;
  flt->SetFunctor(func);
  flt->SetInput1(a);
  flt->SetInput2(b);
  flt->GraftOutput(trg);
  flt->Update();
}

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::compute_semi_lagrangean_a()
{
  for(uint i = 0; i < nt; i++)
    {
    a[i]->FillBuffer(Vec(0.0));
    for (uint j = 0; j < 5; j++)
      {
      interp_vimg(v[i], a[i], -0.5, a[i]);
      vimg_scale_in_place(a[i], dt);
      itk::Index<VDim> x;
      x[0] = 63; x[1] = 63;
      }
    }

}

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::integrate_phi_t0()
{
  // Go through and compute phi_t0 for each time point
  for(int m = 0; m < (int) nt; m++)
    if(m==0)
      {
      f[m]->FillBuffer(Vec(0.0));
      }
    else
      {
      interp_vimg(f[m-1], a[m], -1.0, f[m]);
      vimg_subtract_in_place(f[m], a[m]);
      }
}

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::integrate_phi_t1()
{
  for(int m = nt-1; m >= 0; m--)
    {
    if(m == (int) nt-1)
      {
      f[m]->FillBuffer(Vec(0.0));
      }
    else
      {
      interp_vimg(f[m+1], a[m], 1.0, f[m]);
      vimg_add_in_place(f[m], a[m]);
      }
    } 
}

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::field_jacobian_det(VectorImageType *vec, ImageType *out)
{
  typedef itk::DisplacementFieldJacobianDeterminantFilter<
    VectorImageType, TFloat, ImageType> Filter;
  typename Filter::Pointer filter = Filter::New();
  filter->SetInput(vec);
  filter->SetUseImageSpacingOff();
  filter->GraftOutput(out);
  filter->Update();
}

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::image_gradient(ImageType *src, VectorImageType *grad)
{
  // Create a gradient image filter
  typedef itk::GradientImageFilter<ImageType, TFloat, TFloat> Filter;
  typename Filter::Pointer flt = Filter::New();
  flt->SetInput(src);
  flt->SetUseImageSpacingOff();
  flt->SetUseImageDirection(false);
  flt->GraftOutput(grad);
  flt->Update();
}

template <class TFloat, uint VDim>
void
LDDMMData<TFloat, VDim>
::img_smooth(ImageType *src, ImageType *trg, double sigma)
{
  // typedef itk::SmoothingRecursiveGaussianImageFilter<ImageType, ImageType> Filter;
  typedef itk::DiscreteGaussianImageFilter<ImageType, ImageType> Filter;
  typename Filter::Pointer flt = Filter::New();
  flt->SetInput(src);
  // flt->SetSigma(sigma);
  flt->SetVariance(sigma * sigma);
  flt->GraftOutput(trg);
  flt->Update();
}

template <class TFloat, uint VDim>
void
LDDMMData<TFloat, VDim>
::vimg_smooth(VectorImageType *src, VectorImageType *trg, double sigma)
{
  /// typedef itk::SmoothingRecursiveGaussianImageFilter<VectorImageType, VectorImageType> Filter;
  typedef itk::DiscreteGaussianImageFilter<VectorImageType, VectorImageType> Filter;
  typename Filter::Pointer flt = Filter::New();
  flt->SetInput(src);
  // flt->SetSigma(sigma);
  // flt->GraftOutput(trg);
  flt->SetVariance(sigma * sigma);
  flt->Update();

  vimg_write(src, "preshitme.nii.gz");
  vimg_write(flt->GetOutput(), "shitme.nii.gz");
}




template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::img_read(const char *fn, ImagePointer &trg)
{
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(fn);
  reader->Update();
  trg = reader->GetOutput();
}

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::img_write(ImageType *src, const char *fn)
{
  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(fn);
  writer->SetInput(src);
  writer->Update();
}

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::vimg_read(const char *fn, VectorImagePointer &trg)
{
  typedef itk::ImageFileReader<VectorImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(fn);
  reader->Update();
  trg = reader->GetOutput();
}

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::vimg_write(VectorImageType *src, const char *fn)
{
  // Cast to vector image type
  typedef itk::VectorImage<TFloat, VDim> OutputImageType;
  typename OutputImageType::Pointer output = OutputImageType::New();
  output->CopyInformation(src);
  output->SetRegions(src->GetBufferedRegion());
  output->SetNumberOfComponentsPerPixel(VDim);

  // Override the data pointer
  output->GetPixelContainer()->SetImportPointer(
    (TFloat *) src->GetBufferPointer(), 
    VDim * src->GetPixelContainer()->Size(), false);

  // Write the vector data
  typedef itk::ImageFileWriter<OutputImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(fn);
  writer->SetInput(output);
  writer->Update();

}

template <class TFloat, uint VDim>
void 
LDDMMData<TFloat, VDim>
::vfield_read(uint nt, const char *fnpat, VelocityField &v)
{
  v.clear();
  for(uint i = 0; i < nt; i++)
    {
    char fname[1024];
    sprintf(fname, fnpat, i);
    VectorImagePointer vt;
    vimg_read(fname, vt);
    v.push_back(vt);
    }
}

/* =============================== */

template <class TFloat, uint VDim>
LDDMMFFTInterface<TFloat, VDim>
::LDDMMFFTInterface(ImageType *ref)
{
  // Work out the data dimensions (large enough, and multiple of four)
  m_Size = ref->GetBufferedRegion().GetSize();
  m_Alloc = m_Size;
  m_Alloc[VDim-1] = 2 * (m_Size[VDim-1] / 2 + 1);

  // Size for calling the plan routines
  int n[VDim];

  // Get the data dimensions
  m_AllocSize = 1; m_DataSize = 1;
  for(uint i = 0; i < VDim; i++)
    {
    m_AllocSize *= m_Alloc[i];
    m_DataSize *= m_Size[i];
    n[i] = m_Size[i];
    }

  // Allocate the complex data (real data is packed in the complex data)
  m_Data = (double *) fftw_malloc(m_AllocSize * sizeof(double));

  // Create plans for forward and inverse transforms
  m_Plan = fftw_plan_dft_r2c(VDim, n, m_Data, (fftw_complex *) m_Data, FFTW_MEASURE);
  m_InvPlan = fftw_plan_dft_c2r(VDim, n, (fftw_complex *) m_Data, m_Data, FFTW_MEASURE);
}

template <class TFloat, uint VDim>
void
LDDMMFFTInterface<TFloat, VDim>
::convolution_fft(
  VectorImageType *img, ImageType *kernel_ft, bool inv_kernel,
  VectorImageType *out)
{
  // Pack the data into m_Data. This requires us to skip a few bytes at
  // the end of each row of data
  uint nskip = m_Alloc[VDim-1] - m_Size[VDim-1];
  uint ncopy = m_Size[VDim-1];
  uint nout = m_Alloc[VDim-1] / 2;
  uint noutskip = kernel_ft->GetBufferedRegion().GetSize()[VDim-1] - nout;
  uint nstrides = m_AllocSize / m_Alloc[VDim-1];

  for(uint d = 0; d < VDim; d++)
    {
    const Vec *src = img->GetBufferPointer();
    double *dst = m_Data;

    // Funky loop
    for(double *end = dst + m_AllocSize; dst < end; dst+=nskip)
      for(double *rowend = dst + ncopy; dst < rowend; dst++, src++)
        *dst = (double) (*src)[d];

    // Execute the plan
    fftw_execute(m_Plan);

    // Multiply or divide the complex values in m_Data by the kernel array
    fftw_complex *c = (fftw_complex *) m_Data;

    // Scaling factor for producing final result
    double scale = 1.0 / m_DataSize;

    /*
    // Precision weirdness (results differ from MATLAB fft in 6th, 7th decimal digit)
    uint tp = (m_Alloc[VDim-1] / 2) * 6 + 8;
    printf("Before scaling, value at %d is %12.12lf, %12.12lf\n",
      tp, c[tp][0], c[tp][1]);
    printf("Kernel value at %d is %12.12lf\n", 128*6+8, kp[128*6+8]);
    */

    // Another such loop
    TFloat *kp = kernel_ft->GetBufferPointer();
    if(inv_kernel)
      {
      for(uint i = 0; i < nstrides; i++)
        {
        for(uint j = 0; j < nout; j++)
          {
          (*c)[0] /= (*kp);
          (*c)[1] /= (*kp);
          c++; kp++;
          }
        kp += noutskip;
        }
      }
    else
      {
      for(uint i = 0; i < nstrides; i++)
        {
        for(uint j = 0; j < nout; j++)
          {
          (*c)[0] *= (*kp);
          (*c)[1] *= (*kp);
          c++; kp++;
          }
        kp += noutskip;
        }
      }

    /*
    fftw_complex *ctest = (fftw_complex *) m_Data;
    printf("After scaling, value at %d is %12.12lf, %12.12lf\n",
      tp, ctest[tp][0], ctest[tp][1]);
    */

    // Inverse transform
    fftw_execute(m_InvPlan);

    // Copy the results to the output image
    const double *osrc = m_Data;
    Vec *odst = out->GetBufferPointer();
    for(uint i = 0; i < nstrides; i++, osrc+=nskip)
      for(uint j = 0; j < ncopy; j++, odst++, osrc++)
        (*odst)[d] = (TFloat) ((*osrc) * scale);

    /*
    odst = out->GetBufferPointer();
    printf("Result %12.12lf\n",  odst[128*6+8][0]);
    */
    }

}

template <class TFloat, uint VDim>
LDDMMFFTInterface<TFloat, VDim>
::~LDDMMFFTInterface()
{
  fftw_destroy_plan(m_Plan);
  fftw_destroy_plan(m_InvPlan);
  fftw_free(m_Data);
}

template <class TFloat, uint VDim>
LDDMMImageMatchingObjective<TFloat, VDim>
::LDDMMImageMatchingObjective(LDDMM &p)
  : fft(p.fix)
{
  // Allocate intermediate datasets
  LDDMM::alloc_img(Jt0, p.fix);
  LDDMM::alloc_img(Jt1, p.fix);
  LDDMM::alloc_img(DetPhit1, p.fix);
  LDDMM::alloc_vimg(GradJt0, p.fix);
}

template <class TFloat, uint VDim>
TFloat
LDDMMImageMatchingObjective<TFloat, VDim>
::compute_objective_and_gradient(LDDMM &p)
{
  // Compute the regularization energy of v. We can use a[0] for temporary storage
  // e_field = lddmm_vector_field_dot_product(vx, vy, vx, vy, p);
  TFloat e_field = 0.0;
  for(uint m = 0; m < p.nt; m++)
    {
    // a[0] = Lv[m] .* v[m]
    fft.convolution_fft(p.v[m], p.f_kernel_sq, false, p.a[0]);

    // We're sticking the inner product in Jt0
    LDDMM::vimg_euclidean_inner_product(Jt0, p.a[0], p.v[m]); 
    e_field += LDDMM::img_voxel_sum(Jt0) / p.nt;
    }

  // Compute the 'a' array (for semilagrangean scheme)
  p.compute_semi_lagrangean_a();

  // Go through and compute phi_t1 for each time point
  p.integrate_phi_t1();

  // Compute the update for v at each time step
  for(uint m = 0; m < p.nt; m++)
    {
    // Currently, f[m] holds phi_t1[m]. Use it for whatever we need
    // and then replace with phi_t0[m]

    // TODO: for ft00 and ft11, don't waste time on interpolation

    // Jt1 = lddmm_warp_scalar_field(p.I1, ft1x(:,:,it), ft1y(:,:,it), p);
    LDDMM::interp_img(p.mov, p.f[m], Jt1); 

    // detjac_phi_t1 = lddmm_jacobian_determinant(ft1x(:,:,it), ft1y(:,:,it), p);
    LDDMM::field_jacobian_det(p.f[m], DetPhit1);

    // Place phi_t0 into the f array
    if(m==0)
      {
      p.f[m]->FillBuffer(typename LDDMM::Vec(0.0));
      }
    else
      {
      LDDMM::interp_vimg(p.f[m-1], p.a[m], -1.0, p.f[m]);
      LDDMM::vimg_subtract_in_place(p.f[m], p.a[m]);
      }

    // Jt0 = lddmm_warp_scalar_field(p.I0, ft0x(:,:,it), ft0y(:,:,it), p); 
    LDDMM::interp_img(p.fix, p.f[m], Jt0); 

    // [grad_Jt0_x grad_Jt0_y] = gradient(Jt0);
    LDDMM::image_gradient(Jt0, GradJt0);

    // pde_rhs_x = detjac_phi_t1 .* (Jt0 - Jt1) .* grad_Jt0_x; 
    // pde_rhs_y = detjac_phi_t1 .* (Jt0 - Jt1) .* grad_Jt0_y; 

    // Here we do some small tricks. We want to retain Jt0 because it's the warped
    // template image, and we want to retain the difference Jt0-Jt1 = (J0-I1) for
    // calculating the objective at the end. 
    LDDMM::img_subtract_in_place(Jt1, Jt0);           // 'Jt1' stores Jt1 - Jt0 
    LDDMM::img_multiply_in_place(DetPhit1, Jt1);      // 'DetPhit1' stores (det Phi_t1)(Jt1-Jt0)
    LDDMM::vimg_multiply_in_place(GradJt0, DetPhit1); // 'GradJt0' stores  GradJt0 * (det Phi_t1)(Jt1-Jt0)

    // Solve PDE via FFT convolution
    // pde_soln_x = ifft2(fft2(pde_rhs_x) ./ p.f_kernel_sq,'symmetric');
    // pde_soln_y = ifft2(fft2(pde_rhs_y) ./ p.f_kernel_sq,'symmetric');
    fft.convolution_fft(GradJt0, p.f_kernel_sq, true, GradJt0); // 'GradJt0' stores K[ GradJt0 * (det Phi_t1)(Jt1-Jt0) ]

    // dedvx(:,:,it) = dedvx(:,:,it) - 2 * pde_soln_x / p.sigma^2;
    // dedvy(:,:,it) = dedvy(:,:,it) - 2 * pde_soln_y / p.sigma^2;        

    // Store the update in a[m]
    LDDMM::vimg_scale_in_place(GradJt0, 1.0 / p.sigma_sq); // 'GradJt0' stores 1 / sigma^2 K[ GradJt0 * (det Phi_t1)(Jt1-Jt0) ]
    LDDMM::vimg_add_in_place(GradJt0, p.v[m]); // 'GradJt0' stores v + 1 / sigma^2 K[ GradJt0 * (det Phi_t1)(Jt1-Jt0) ]
    LDDMM::vimg_scale(GradJt0, 2.0, p.a[m]); // p.a[m] holds 2 v + 2 / sigma^2 K[ GradJt0 * (det Phi_t1)(Jt1-Jt0) ]
    }

  // Ok, Jt1 currently contains (Jt1-Jt0), we just need to square it.
  TFloat e_image = LDDMM::img_euclidean_norm_sq(Jt1) / p.sigma_sq;

  // Return the energy
  printf("  Energy components: %lf, %lf\n", e_field, e_image);
  return e_field + e_image;
}


/*
template class LDDMMData<float, 2>;
template class LDDMMData<float, 3>;
template class LDDMMFFTInterface<float, 2>;
template class LDDMMFFTInterface<float, 3>;
*/

template class LDDMMData<double, 2>;
template class LDDMMData<double, 3>;
template class LDDMMData<double, 4>;

template class LDDMMFFTInterface<myreal, 2>;
template class LDDMMFFTInterface<myreal, 3>;
template class LDDMMImageMatchingObjective<myreal, 2>;
template class LDDMMImageMatchingObjective<myreal, 3>;

