#ifndef _LDDMM_DATA_H_
#define _LDDMM_DATA_H_

#include <lddmm_common.h>
#include <itkNumericTraits.h>
#include <itkImage.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkCovariantVector.h>
#include <vnl/vnl_math.h>
#include <vector>

#include <fftw3.h>

template<class TFloat, uint VDim>
class LDDMMData
{
public:
  // Image data
  typedef itk::ImageBase<VDim> ImageBaseType;
  typedef itk::Image<TFloat, VDim> ImageType;
  typedef typename ImageType::Pointer ImagePointer;
  typedef itk::ImageRegionIteratorWithIndex<ImageType> ImageIterator;

  // Vector fields, etc
  typedef itk::CovariantVector<TFloat, VDim> Vec;
  typedef itk::Image<Vec, VDim> VectorImageType;
  typedef typename VectorImageType::Pointer VectorImagePointer;
  typedef std::vector<VectorImagePointer> VelocityField;

  // Regions, etc
  typedef itk::ImageRegion<VDim> RegionType;

  // Pointers to the fixed and moving images
  ImagePointer fix, mov;

  // Fourier space kernels
  ImagePointer f_kernel, f_kernel_sq;

  // Velocity field pointers (v, phi, a used for semi-lagrange scheme)
  VelocityField v, f, a;

  // Region for the velocity fields
  RegionType r;

  // Parameters
  double alpha, sigma, gamma, dt, sigma_sq;
  
  // Dimensions
  uint n[VDim];

  // Number of timesteps, number of voxels
  uint nt, nv;

  // Allocate a velocity field
  static void alloc_vf(VelocityField &vf, uint nt, ImageBaseType *ref);
  static void alloc_img(ImagePointer &img, ImageBaseType *ref);
  static void alloc_vimg(VectorImagePointer &vimg, ImageBaseType *ref);

  // Initialize LDDMM data 
  static void init(LDDMMData<TFloat, VDim> &, 
    ImageType *fix, ImageType *mov, 
    uint nt, double alpha, double gamma, double sigma);

  // Apply deformation to data
  static void interp_vimg(
    VectorImageType *data, VectorImageType *field, 
    TFloat def_scale, VectorImageType *out);

  // Apply deformation to data
  static void interp_img(ImageType *data, VectorImageType *field, ImageType *out);

  // Take Jacobian of deformation field
  static void field_jacobian_det(VectorImageType *vec, ImageType *out);

  // Smooth an image in-place
  static void img_smooth(ImageType *src, ImageType *out, double sigma);
  static void vimg_smooth(VectorImageType *src, VectorImageType *out, double sigma);

  // Take gradient of an image
  static void image_gradient(ImageType *src, VectorImageType *grad);

  // Basic math
  static void vimg_add_in_place(VectorImagePointer &trg, VectorImageType *a);
  static void vimg_subtract_in_place(VectorImagePointer &trg, VectorImageType *a);
  static void vimg_scale_in_place(VectorImagePointer &trg, TFloat s);

  // compute trg = trg + s * a
  static void vimg_add_scaled_in_place(VectorImagePointer &trg, VectorImageType *a, TFloat s);

  static void vimg_scale(VectorImageType *src, TFloat s, VectorImagePointer &trg);
  static void vimg_multiply_in_place(VectorImagePointer &trg, ImageType *s);
  static void vimg_euclidean_inner_product(ImagePointer &trg, VectorImageType *a, VectorImageType *b);
  static TFloat vimg_euclidean_norm_sq(VectorImageType *trg);

  // Scalar math
  static void img_subtract_in_place(ImagePointer &trg, ImageType *a);
  static void img_multiply_in_place(ImagePointer &trg, ImageType *a);
  static TFloat img_euclidean_norm_sq(ImageType *trg);
  static TFloat img_voxel_sum(ImageType *trg);
  static void img_min_max(ImageType *src, TFloat &out_min, TFloat &out_max);
  
  // Some IO methods
  static void img_read(const char *fn, ImagePointer &trg);
  static void img_write(ImageType *src, const char *fn);
  static void vimg_read(const char *fn, VectorImagePointer &trg);
  static void vimg_write(VectorImageType *src, const char *fn);

  static void vfield_read(uint nt, const char *fnpat, VelocityField &v);

  // Compute a array from v
  void compute_semi_lagrangean_a();

  // Integrate forward tranform (phi_0_t)
  void integrate_phi_t0();
  void integrate_phi_t1();

protected:

  // A vector image for in-place interpolation operations
  VectorImagePointer vtmp;

};

template <class TFloat, uint VDim>
class LDDMMFFTInterface
{
public:
  typedef typename LDDMMData<TFloat, VDim>::ImageType ImageType;
  typedef typename LDDMMData<TFloat, VDim>::VectorImageType VectorImageType;
  typedef typename LDDMMData<TFloat, VDim>::Vec Vec;

  LDDMMFFTInterface(ImageType *ref);
  ~LDDMMFFTInterface();

  void convolution_fft(
    VectorImageType *img, ImageType *kernel_ft, bool inv_kernel, 
    VectorImageType *out);

private:

  // Size of the input array and allocated array (bigger, for in-place math)
  itk::Size<VDim> m_Size, m_Alloc;
  uint m_AllocSize, m_DataSize;

  // In-place data array
  double *m_Data;

  // FFT plan
  fftw_plan m_Plan, m_InvPlan;

};

// Class for iteratively computing the objective function
template<class TFloat, uint VDim>
class LDDMMImageMatchingObjective 
{
public:
  typedef LDDMMData<TFloat, VDim> LDDMM;
  typedef LDDMMFFTInterface<TFloat, VDim> FFT;

  LDDMMImageMatchingObjective(LDDMM &p);
  TFloat compute_objective_and_gradient(LDDMM &p);

  typename LDDMM::ImagePointer Jt0, Jt1, DetPhit1;
  typename LDDMM::VectorImagePointer GradJt0;
  FFT fft;
};




#endif
