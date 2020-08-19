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

#include <iostream>
#include <sstream>
#include <cstdio>
#include <vector>
#include <string>
#include <algorithm>

#include "lddmm_common.h"
#include "lddmm_data.h"

#include <itkImageFileReader.h>
#include <itkAffineTransform.h>
#include <itkTransformFactory.h>
#include <itkTimeProbe.h>
#include <itkImageFileWriter.h>
#include <itkCastImageFilter.h>

#include "MultiImageRegistrationHelper.h"
#include "FastWarpCompositeImageFilter.h"
#include "MultiComponentImageMetricBase.h"

#include <vnl/algo/vnl_powell.h>
#include <vnl/algo/vnl_svd.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vnl/vnl_trace.h>
#include <vnl/vnl_numeric_traits.h>



// Helper function to get the RAS coordinate of the center of 
// an image
template <unsigned int VDim>
vnl_vector<double>
GetImageCenterinNiftiSpace(itk::ImageBase<VDim> *image)
{
  itk::ImageRegion<VDim> r = image->GetBufferedRegion();
  itk::ContinuousIndex<double, VDim> idx;
  itk::Point<double, VDim> ctr;
  for(int d = 0; d < VDim; d++)
    idx[d] = r.GetIndex()[d] + r.GetSize()[d] * 0.5;
  image->TransformContinuousIndexToPhysicalPoint(idx, ctr);

  // Map to RAS (flip first two coordinates)
  for(int d = 0; d < 2 && d < VDim; d++)
    ctr[d] = -ctr[d];

  return ctr.GetVnlVector();
}




#include "itkTransformFileReader.h"

template <unsigned int VDim, typename TReal>
vnl_matrix<double>
GreedyApproach<VDim, TReal>
::ReadAffineMatrixViaCache(const TransformSpec &ts)
{
  // Physical (RAS) space transform matrix
  vnl_matrix<double> Qp(VDim+1, VDim+1);  Qp.set_identity();

  // An ITK-style transform - forced to floating point here
  typedef itk::MatrixOffsetTransformBase<double, VDim, VDim> TransformType;
  typename TransformType::Pointer itk_tran;

  // See if a transform is already stored in the cache
  typename ImageCache::const_iterator itCache = m_ImageCache.find(ts.filename);
  if(itCache != m_ImageCache.end())
    {
    TransformType *cached = dynamic_cast<TransformType *>(itCache->second.target);
    if(!cached)
      throw GreedyException("Cached transform %s cannot be cast to type %s",
                            ts.filename.c_str(), typeid(TransformType).name());

    itk_tran = cached;
    }
  else
    {
    // Open the file and read the first line
    std::ifstream fin(ts.filename.c_str());
    std::string header_line, itk_header = "#Insight Transform File";
    std::getline(fin, header_line);

    if(header_line.substr(0, itk_header.size()) == itk_header)
      {
      fin.close();
      try
        {
        // First we try to load the transform using ITK code
        // This code is from c3d_affine_tool
        typedef itk::AffineTransform<double, VDim> AffTran;
        itk::TransformFactory<TransformType>::RegisterTransform();
        itk::TransformFactory<AffTran>::RegisterTransform();

        itk::TransformFileReader::Pointer fltReader = itk::TransformFileReader::New();
        fltReader->SetFileName(ts.filename.c_str());
        fltReader->Update();

        itk::TransformBase *base = fltReader->GetTransformList()->front();
        itk_tran = dynamic_cast<TransformType *>(base);
        }
      catch(...)
        {
        throw GreedyException("Unable to read ITK transform file %s", ts.filename.c_str());
        }
      }
    else
      {
      // Try reading C3D matrix format
      fin.seekg(0);
      for(size_t i = 0; i < VDim+1; i++)
        for(size_t j = 0; j < VDim+1; j++)
          if(fin.good())
            {
            fin >> Qp[i][j];
            }
      fin.close();
      }
    }

  // At this point we might have read the RAS matrix directly, or an ITK transform
  // if the latter, extract the RAS matrix
  if(itk_tran.IsNotNull())
    {
    for(size_t r = 0; r < VDim; r++)
      {
      for(size_t c = 0; c < VDim; c++)
        {
        Qp(r,c) = itk_tran->GetMatrix()(r,c);
        }
      Qp(r,VDim) = itk_tran->GetOffset()[r];
      }

    // RAS - LPI nonsense
    if(VDim == 3)
      {
      Qp(2,0) *= -1; Qp(2,1) *= -1;
      Qp(0,2) *= -1; Qp(1,2) *= -1;
      Qp(0,3) *= -1; Qp(1,3) *= -1;
      }
    }

  // Compute the exponent
  if(ts.exponent == 1.0)
    {
    return Qp;
    }
  else if(ts.exponent == -1.0)
    {
    return vnl_matrix_inverse<double>(Qp);
    }
  else
    {
    throw GreedyException("Transform exponent values of +1 and -1 are the only ones currently supported");
    }

  return Qp;
}




template <unsigned int VDim, typename TReal>
void
GreedyApproach<VDim, TReal>
::WriteAffineMatrixViaCache(
    const std::string &filename, const vnl_matrix<double> &Qp)
{
  // An ITK-style transform - forced to double point here
  typedef itk::MatrixOffsetTransformBase<double, VDim, VDim> TransformType;

  // See if a transform is already stored in the cache
  typename ImageCache::const_iterator itCache = m_ImageCache.find(filename);
  if(itCache != m_ImageCache.end())
    {
    TransformType *cached = dynamic_cast<TransformType *>(itCache->second.target);
    if(!cached)
      throw GreedyException("Cached transform %s cannot be cast to type %s",
                            filename.c_str(), typeid(TransformType).name());

    // RAS - LPI nonsense
    vnl_matrix<double> Q = Qp;
    if(VDim == 3)
      {
      Q(2,0) *= -1; Q(2,1) *= -1;
      Q(0,2) *= -1; Q(1,2) *= -1;
      Q(0,3) *= -1; Q(1,3) *= -1;
      }

    typename TransformType::MatrixType matrix;
    typename TransformType::OffsetType offset;

    // We have found the output transform and can use it for assignment
    for(size_t r = 0; r < VDim; r++)
      {
      for(size_t c = 0; c < VDim; c++)
        {
        matrix(r, c) = Q(r, c);
        }
      offset[r] = Q(r, VDim);
      }

    cached->SetMatrix(matrix);
    cached->SetOffset(offset);
    }

  // Write to actual file
  if(itCache == m_ImageCache.end() || itCache->second.force_write)
    {
    std::ofstream matrixFile;
    matrixFile.open(filename.c_str());
    matrixFile << Qp;
    matrixFile.close();
    }
}

template<unsigned int VDim, typename TReal>
vnl_matrix<double>
GreedyApproach<VDim, TReal>
::ReadAffineMatrix(const TransformSpec &ts)
{
  GreedyApproach<VDim, TReal> api;
  return api.ReadAffineMatrixViaCache(ts);
}

template<unsigned int VDim, typename TReal>
void
GreedyApproach<VDim, TReal>
::WriteAffineMatrix(const std::string &filename, const vnl_matrix<double> &Qp)
{
  GreedyApproach<VDim, TReal> api;
  api.WriteAffineMatrixViaCache(filename, Qp);
}

template <unsigned int VDim, typename TReal>
template <class TImage>
itk::SmartPointer<TImage>
GreedyApproach<VDim, TReal>
::ReadImageViaCache(const std::string &filename,
                    itk::ImageIOBase::IOComponentType *comp_type)
{
  // Check the cache for the presence of the image
  typename ImageCache::const_iterator it = m_ImageCache.find(filename);
  if(it != m_ImageCache.end())
    {
    itk::Object *cached_object = it->second.target;
    TImage *image = dynamic_cast<TImage *>(cached_object);
    if(!image)
      throw GreedyException("Cached image %s cannot be cast to type %s",
                            filename.c_str(), typeid(TImage).name());
    itk::SmartPointer<TImage> pointer = image;

    // The component type is unknown here
    if(comp_type)
      *comp_type = itk::ImageIOBase::UNKNOWNCOMPONENTTYPE;

    return pointer;
    }

  // Read the image using ITK reader
  typedef itk::ImageFileReader<TImage> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(filename.c_str());
  reader->Update();

  // Store the component type if requested
  if(comp_type)
    *comp_type = reader->GetImageIO()->GetComponentType();

  itk::SmartPointer<TImage> pointer = reader->GetOutput();
  return pointer;
}

template <unsigned int VDim, typename TReal>
typename GreedyApproach<VDim, TReal>::ImageBaseType::Pointer
GreedyApproach<VDim, TReal>
::ReadImageBaseViaCache(const std::string &filename)
{
  // Check the cache for the presence of the image
  typename ImageCache::const_iterator it = m_ImageCache.find(filename);
  if(it != m_ImageCache.end())
    {
    ImageBaseType *image_base = dynamic_cast<ImageBaseType *>(it->second.target);
    if(!image_base)
      throw GreedyException("Cached image %s cannot be cast to type %s",
                            filename.c_str(), typeid(ImageBaseType).name());
    typename ImageBaseType::Pointer pointer = image_base;
    return pointer;
    }

  // Read the image using ITK reader
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(filename.c_str());
  reader->Update();

  typename ImageBaseType::Pointer pointer = reader->GetOutput();
  return pointer;
}


template <unsigned int VDim, typename TReal>
template <class TImage>
void
GreedyApproach<VDim, TReal>
::WriteImageViaCache(TImage *img, const std::string &filename, typename LDDMMType::IOComponentType comp)
{
  typename ImageCache::const_iterator it = m_ImageCache.find(filename);
  if(it != m_ImageCache.end())
    {
    // The image was found in the cache. Make sure it is an image pointer
    typename LDDMMType::ImageBaseType *cached =
        dynamic_cast<typename LDDMMType::ImageBaseType *>(it->second.target);

    if(!cached)
      throw GreedyException("Cached image %s cannot be cast to ImageBase",
                            filename.c_str(), typeid(TImage).name());

    // Try the automatic cast
    bool cast_rc = false;

    // This is a little dumb, but the vimg_write uses some special code that
    // we need to account for
    if(dynamic_cast<VectorImageType *>(img))
      cast_rc =LDDMMType::vimg_auto_cast(dynamic_cast<VectorImageType *>(img), cached);
    else if(dynamic_cast<ImageType *>(img))
      cast_rc =LDDMMType::img_auto_cast(dynamic_cast<ImageType *>(img), cached);
    else if(dynamic_cast<CompositeImageType *>(img))
      cast_rc =LDDMMType::cimg_auto_cast(dynamic_cast<CompositeImageType *>(img), cached);
    else
      {
      // Some other type (e.g., LabelImage). Instead of doing an auto_cast, we require the
      // cached object to be of the same type as the image
      TImage *cached_typed = dynamic_cast<TImage *>(cached);
      if(cached_typed)
        {
        typedef itk::CastImageFilter<TImage, TImage> CopyFilterType;
        typename CopyFilterType::Pointer copier = CopyFilterType::New();
        copier->SetInput(img);
        copier->GraftOutput(cached_typed);
        copier->Update();
        }
      else throw GreedyException("Cached image %s cannot be cast to type %s",
                                 filename.c_str(), typeid(TImage).name());
      }


    // If cast failed, throw exception
    if(!cast_rc)
      throw GreedyException("Image to save %s could not cast to any known type", filename.c_str());
    }

  if(it == m_ImageCache.end() || it->second.force_write)
    {
    // This is a little dumb, but the vimg_write uses some special code that
    // we need to account for
    if(dynamic_cast<VectorImageType *>(img))
      LDDMMType::vimg_write(dynamic_cast<VectorImageType *>(img), filename.c_str(), comp);
    else if(dynamic_cast<ImageType *>(img))
      LDDMMType::img_write(dynamic_cast<ImageType *>(img), filename.c_str(), comp);
    else if(dynamic_cast<CompositeImageType *>(img))
      LDDMMType::cimg_write(dynamic_cast<CompositeImageType *>(img), filename.c_str(), comp);
    else
      {
      // Some other type (e.g., LabelImage). We use the image writer and ignore the comp
      typedef itk::ImageFileWriter<TImage> WriterType;
      typename WriterType::Pointer writer = WriterType::New();
      writer->SetFileName(filename.c_str());
      writer->SetUseCompression(true);
      writer->SetInput(img);
      writer->Update();
      }
    }
}



#include <itkBinaryErodeImageFilter.h>

template <unsigned int VDim, typename TReal>
void GreedyApproach<VDim, TReal>
::ReadImages(GreedyParameters &param, OFHelperType &ofhelper)
{
  // If the parameters include a sequence of transforms, apply it first
  VectorImagePointer moving_pre_warp;

  // Keep a pointer to the fixed image space
  typename OFHelperType::ImageBaseType *ref_space = NULL;

  // Read the input images and stick them into an image array
  for(int i = 0; i < param.inputs.size(); i++)
    {
    // Read fixed and moving images
    CompositeImagePointer imgFix = ReadImageViaCache<CompositeImageType>(param.inputs[i].fixed);
    CompositeImagePointer imgMov = ReadImageViaCache<CompositeImageType>(param.inputs[i].moving);

    // Store the fixed image and/or check it
    if(ref_space == NULL)
      ref_space = imgFix;

    // Read the pre-warps (only once)
    if(param.moving_pre_transforms.size() && moving_pre_warp.IsNull())
      {
      ReadTransformChain(param.moving_pre_transforms, imgFix, moving_pre_warp);
      }

    if(moving_pre_warp.IsNotNull())
      {
      // Create an image to store the warp
      CompositeImagePointer warped_moving =
          LDDMMType::new_cimg(imgFix, imgMov->GetNumberOfComponentsPerPixel());

      // Interpolate the moving image using the transform chain
      LDDMMType::interp_cimg(imgMov, moving_pre_warp, warped_moving, false, true);

      // Add the image pair to the helper
      ofhelper.AddImagePair(imgFix, warped_moving, param.inputs[i].weight);
      }
    else
      {
      // Add to the helper object
      ofhelper.AddImagePair(imgFix, imgMov, param.inputs[i].weight);
      }
    }

  // Read the fixed-space mask
  if(param.gradient_mask.size())
    {
    typedef typename OFHelperType::FloatImageType MaskType;
    typename MaskType::Pointer imgMask =
        ReadImageViaCache<MaskType>(param.gradient_mask);

    ofhelper.SetGradientMask(imgMask);
    }

  if(param.gradient_mask_trim_radius.size() == VDim)
    {
    if(param.gradient_mask.size())
      throw GreedyException("Cannot specify both gradient mask and gradient mask trim radius");

    ofhelper.SetGradientMaskTrimRadius(param.gradient_mask_trim_radius);
    }

  // Read the moving-space mask
  if(param.moving_mask.size())
    {
    typedef typename OFHelperType::FloatImageType MaskType;
    typename MaskType::Pointer imgMovMask =
        ReadImageViaCache<MaskType>(param.moving_mask);

    if(moving_pre_warp.IsNotNull())
      {
      // Create an image to store the warp
      typename MaskType::Pointer warped_moving_mask = LDDMMType::new_img(moving_pre_warp);

      // Interpolate the moving image using the transform chain
      LDDMMType::interp_img(imgMovMask, moving_pre_warp, warped_moving_mask, false, true);

      // Add the warped mask to the helper
      ofhelper.SetMovingMask(warped_moving_mask);
      }
    else
      {
      // Add the mask to the helper object
      ofhelper.SetMovingMask(imgMovMask);
      }
    }

  // Generate the optimized composite images. For the NCC metric, we add random noise to
  // the composite images, specified in units of the interquartile intensity range.
  double noise = (param.metric == GreedyParameters::NCC) ? param.ncc_noise_factor : 0.0;

  // Build the composite images
  ofhelper.BuildCompositeImages(noise);

  // If the metric is NCC, then also apply special processing to the gradient masks
  if(param.metric == GreedyParameters::NCC)
    ofhelper.DilateCompositeGradientMasksForNCC(array_caster<VDim>::to_itkSize(param.metric_radius));
}

#include <vnl/algo/vnl_lbfgs.h>

template <unsigned int VDim, typename TReal>
vnl_matrix<double>
GreedyApproach<VDim, TReal>
::MapAffineToPhysicalRASSpace(
    OFHelperType &of_helper, int level,
    LinearTransformType *tran)
{
  // Map the transform to NIFTI units
  vnl_matrix<double> T_fix, T_mov, Q, A;
  vnl_vector<double> s_fix, s_mov, p, b;

  GetVoxelSpaceToNiftiSpaceTransform(of_helper.GetReferenceSpace(level), T_fix, s_fix);
  GetVoxelSpaceToNiftiSpaceTransform(of_helper.GetMovingReferenceSpace(level), T_mov, s_mov);

  itk_matrix_to_vnl_matrix(tran->GetMatrix(), A);
  itk_vector_to_vnl_vector(tran->GetOffset(), b);

  Q = T_mov * A * vnl_matrix_inverse<double>(T_fix);
  p = T_mov * b + s_mov - Q * s_fix;

  vnl_matrix<double> Qp(VDim+1, VDim+1);
  Qp.set_identity();
  for(int i = 0; i < VDim; i++)
    {
    Qp(i, VDim) = p(i);
    for(int j = 0; j < VDim; j++)
      Qp(i,j) = Q(i,j);
    }

  return Qp;
}

template <unsigned int VDim, typename TReal>
void
GreedyApproach<VDim, TReal>
::MapPhysicalRASSpaceToAffine(
    OFHelperType &of_helper, int level,
    vnl_matrix<double> &Qp,
    LinearTransformType *tran)
{
  // Map the transform to NIFTI units
  vnl_matrix<double> T_fix, T_mov, Q(VDim, VDim), A;
  vnl_vector<double> s_fix, s_mov, p(VDim), b;

  GetVoxelSpaceToNiftiSpaceTransform(of_helper.GetReferenceSpace(level), T_fix, s_fix);
  GetVoxelSpaceToNiftiSpaceTransform(of_helper.GetMovingReferenceSpace(level), T_mov, s_mov);

  for(int i = 0; i < VDim; i++)
    {
    p(i) = Qp(i, VDim);
    for(int j = 0; j < VDim; j++)
      Q(i,j) = Qp(i,j);
    }

  // A = vnl_matrix_inverse<double>(T_mov) * (Q * T_fix);
  // b = vnl_matrix_inverse<double>(T_mov) * (p - s_mov + Q * s_fix);
  A=vnl_svd<double>(T_mov).solve(Q * T_fix);
  b=vnl_svd<double>(T_mov).solve(p - s_mov + Q * s_fix);

  typename LinearTransformType::MatrixType tran_A;
  typename LinearTransformType::OffsetType tran_b;

  vnl_matrix_to_itk_matrix(A, tran_A);
  vnl_vector_to_itk_vector(b, tran_b);

  tran->SetMatrix(tran_A);
  tran->SetOffset(tran_b);
}

template <unsigned int VDim, typename TReal>
void
GreedyApproach<VDim, TReal>
::RecordMetricValue(const MultiComponentMetricReport &metric)
{
  if(m_MetricLog.size())
    m_MetricLog.back().push_back(metric);
}

/**
 * Find a plane of symmetry in an image
 */
/*
template <unsigned int VDim, typename TReal>
vnl_vector<double>
GreedyApproach<VDim, TReal>
::FindSymmetryPlane(ImageType *image, int N, int n_search_pts)
{
  typedef vnl_vector_fixed<double, 3> Vec3;
  typedef vnl_matrix_fixed<double, 3, 3> Mat3;

  // Loop over direction on a sphere, using the Saff & Kuijlaars algorithm
  // https://perswww.kuleuven.be/~u0017946/publications/Papers97/art97a-Saff-Kuijlaars-MI/Saff-Kuijlaars-MathIntel97.pdf
  double phi = 0.0;
  double spiral_const = 3.6 / sqrt(N);
  for(int k = 0; k < n_sphere_pts; k++)
    {
    // Height of the k-th point
    double cos_theta = -1 * (2 * k) / (N - 1);
    double sin_theta = sqrt(1 - cos_theta * cos_theta);

    // Phase of the k-th point
    if(k > 0 && k < N-1)
      phi = fmod(phi_last + spiral_const / sin_theta, vnl_math::pi * 2);
    else
      phi = 0.0;

    // We now have the polar coordinates of the points, get cartesian coordinates
    Vec3 q;
    q[0] = sin_theta * cos(phi);
    q[1] = sin_theta * sin(phi);
    q[2] = cos_theta;

    // Now q * (x,y,z) = 0 defines a plane through the origin. We will test whether the image
    // is symmetric across this plane. We first construct the reflection matrix
    Mat3 R;
    R(0,0) =  1 - q[0] * q[0]; R(0,1) = -2 * q[1] * q[0]; R(0,2) = -2 * q[2] * q[0];
    R(1,0) = -2 * q[0] * q[1]; R(1,1) =  1 - q[1] * q[1]; R(1,2) = -2 * q[2] * q[1];
    R(2,0) = -2 * q[0] * q[2]; R(2,1) = -2 * q[1] * q[2]; R(2,2) =  1 - q[2] * q[2];

    // We must find the reasonable range of intercepts to search for. An intercept is reasonable
    // if the plane cuts the image volume in at least a 80/20 ratio (let's say)


    // This is a test axis of rotation. We will now measure the symmetry of the image across this axis
    // To do so, we will construct a series of flips across this direction

    }
}
*/

/**
 * This method performs initial alignment by first searching for a plane of symmetry
 * in each image, and then finding the transform between the planes of symmetry.
 *
 * The goal is to have an almost sure-fire method for brain alignment, yet generic
 * enough to work for other data as well.
 */
/*
template <unsigned int VDim, typename TReal>
int GreedyApproach<VDim, TReal>
::SymmetrySearch(GreedyParameters &param, int level, OFHelperType *of_helper)
{

}
*/


template <unsigned int VDim, typename TReal>
int GreedyApproach<VDim, TReal>
::RunAffine(GreedyParameters &param)
{
  typedef AbstractAffineCostFunction<VDim, TReal> AbstractAffineCostFunction;
  typedef RigidCostFunction<VDim, TReal> RigidCostFunction;
  typedef ScalingCostFunction<VDim, TReal> ScalingCostFunction;
  typedef PhysicalSpaceAffineCostFunction<VDim, TReal> PhysicalSpaceAffineCostFunction;

  // Create an optical flow helper object
  OFHelperType of_helper;

  // Set the scaling factors for multi-resolution
  of_helper.SetDefaultPyramidFactors(param.iter_per_level.size());

  // Add random sampling jitter for affine stability at voxel edges
  of_helper.SetJitterSigma(param.affine_jitter);

  // Read the image pairs to register - this will also build the composite pyramids
  ReadImages(param, of_helper);

  // Matrix describing current transform in physical space
  vnl_matrix<double> Q_physical;

  // The number of resolution levels
  unsigned nlevels = param.iter_per_level.size();

  // Clear the metric log
  m_MetricLog.clear();

  // Iterate over the resolution levels
  for(unsigned int level = 0; level < nlevels; ++level)
    {
    // Add stage to metric log
    m_MetricLog.push_back(std::vector<MultiComponentMetricReport>());

    // Define the affine cost function
    AbstractAffineCostFunction *pure_acf, *acf;
    if(param.affine_dof == GreedyParameters::DOF_RIGID)
      {
      RigidCostFunction *rigid_acf = new RigidCostFunction(&param, this, level, &of_helper);
      acf = new ScalingCostFunction(
              rigid_acf,
              rigid_acf->GetOptimalParameterScaling(
                of_helper.GetReferenceSpace(level)->GetBufferedRegion().GetSize()));
      pure_acf = rigid_acf;
      }
    else
      {
      //  PureAffineCostFunction *affine_acf = new PureAffineCostFunction(&param, level, &of_helper);
      PhysicalSpaceAffineCostFunction *affine_acf = new PhysicalSpaceAffineCostFunction(&param, this, level, &of_helper);
      acf = new ScalingCostFunction(
              affine_acf,
              affine_acf->GetOptimalParameterScaling(
                of_helper.GetReferenceSpace(level)->GetBufferedRegion().GetSize()));
      pure_acf = affine_acf;
      }

    // Current transform
    typename LinearTransformType::Pointer tLevel = LinearTransformType::New();

    // Set up the initial transform
    if(level == 0)
      {
      // Get the coefficients corresponding to the identity transform in voxel space
      tLevel->SetIdentity();
      vnl_vector<double> xIdent = acf->GetCoefficients(tLevel);

      // Use the provided initial affine as the starting point
      if(param.affine_init_mode == RAS_FILENAME)
        {
        // Read the initial affine transform from a file
        vnl_matrix<double> Qp = this->ReadAffineMatrixViaCache(param.affine_init_transform);

        // Map this to voxel space
        MapPhysicalRASSpaceToAffine(of_helper, level, Qp, tLevel);
        }
      else if(param.affine_init_mode == RAS_IDENTITY)
        {
        // Physical space transform
        vnl_matrix<double> Qp(VDim+1, VDim+1); Qp.set_identity();

        // Map this to voxel space
        MapPhysicalRASSpaceToAffine(of_helper, level, Qp, tLevel);
        }
      else if(param.affine_init_mode == IMG_CENTERS)
        {
        // Find a translation that maps center voxel of fixed image to the center 
        // voxel of the moving image
        vnl_matrix<double> Qp(VDim+1, VDim+1); Qp.set_identity();
        vnl_vector<double> cfix = GetImageCenterinNiftiSpace(of_helper.GetReferenceSpace(level));
        vnl_vector<double> cmov = GetImageCenterinNiftiSpace(of_helper.GetMovingReferenceSpace(level));

        // TODO: I think that setting the matrix above to affine will break the registration
        // if fixed and moving are in different orientations? Or am I crazy?

        // Compute the transform that takes fixed into moving
        for(int d = 0; d < VDim; d++)
          Qp(d, VDim) = cmov[d] - cfix[d];

        // Map this to voxel space
        MapPhysicalRASSpaceToAffine(of_helper, level, Qp, tLevel);
        }

      // Get the new coefficients
      vnl_vector<double> xInit = acf->GetCoefficients(tLevel);

      // If the voxel-space transform is identity, apply a little bit of jitter
      if((xIdent - xInit).inf_norm() < 1e-4)
        {
        // Apply jitter
        vnl_random rndy(12345);
        for (unsigned i = 0; i < xInit.size(); i++)
          xInit[i] += rndy.drand32(-0.4, 0.4);

        // Map back into transform format
        acf->GetTransform(xInit, tLevel);
        }

      // If the uses asks for rigid search, do it!
      if(param.rigid_search.iterations > 0)
        {
        // Random seed. TODO: let user supply seed
        vnl_random randy(12345);

        // For rigid search, we must search in physical space, rather than in voxel space.
        // This is the affine transformation in physical space that corresponds to whatever
        // the current initialization is.
        vnl_matrix<double> Qp = MapAffineToPhysicalRASSpace(of_helper, level, tLevel);

        // Get the center of the fixed image in physical coordinates
        vnl_vector<double> cfix = GetImageCenterinNiftiSpace(of_helper.GetReferenceSpace(level));

        // Create a pure rigid acf
        RigidCostFunction search_fun(&param, this, level, &of_helper);

        // Report the initial best
        double fBest = 0.0;
        vnl_vector<double> xBest = search_fun.GetCoefficients(tLevel);
        search_fun.compute(xBest, &fBest, NULL);
        std::cout << "Rigid search -> Initial best: " << fBest << " " << xBest << std::endl;

        // Loop over random iterations
        for(int i = 0; i < param.rigid_search.iterations; i++)
          {
          // Depending on the search mode, we either apply a small rotation, or any random rotation,
          // or a random rotation and a flip to the input. Whatever rotation we apply, it must
          // be around the center of the fixed coordinate system.
          typename RigidCostFunction::Mat RF;
          if(param.rigid_search.mode == RANDOM_NORMAL_ROTATION)
            {
            // Random angle in radians
            double alpha = randy.normal() * param.rigid_search.sigma_angle * 0.01745329252;
            RF = RigidCostFunction::GetRandomRotation(randy, alpha);
            }
          else if(param.rigid_search.mode == ANY_ROTATION)
            {
            double alpha = randy.drand32(-vnl_math::pi, vnl_math::pi);
            RF = RigidCostFunction::GetRandomRotation(randy, alpha);
            }
          else if(param.rigid_search.mode == ANY_ROTATION_AND_FLIP)
            {
            typename RigidCostFunction::Mat R, F;
            F.set_identity();
            for(unsigned int a = 0; a < VDim; a++)
              F(a,a) = (randy.normal() > 0.0) ? 1.0 : -1.0;
            double alpha = randy.drand32(-vnl_math::pi, vnl_math::pi);
            R = RigidCostFunction::GetRandomRotation(randy, alpha);
            RF = R * F;
            }
          else throw GreedyException("Unknown rotation search mode encountered");

          // Find the offset so that the rotation/flip preserve fixed image center
          typename RigidCostFunction::Vec b_RF = cfix - RF * cfix;

          // Create the physical space matrix corresponding to random search point
          vnl_matrix<double> Qp_rand(VDim+1, VDim+1); Qp_rand.set_identity();
          Qp_rand.update(RF);
          for(unsigned int a = 0; a < VDim; a++)
            Qp_rand(a,VDim) = b_RF[a];

          // Combine the two matrices. The matrix Qp_rand operates in fixed image space so
          // it should be applied first, followed by Qp
          vnl_matrix<double> Qp_search = Qp * Qp_rand;

          // Add the random translation
          for(unsigned int a = 0; a < VDim; a++)
            Qp_search(a,VDim) += randy.normal() * param.rigid_search.sigma_xyz;

          // Convert this physical space transformation into a voxel-space transform
          typename LinearTransformType::Pointer tSearchTry = LinearTransformType::New();
          MapPhysicalRASSpaceToAffine(of_helper, level, Qp_search, tSearchTry);

          // Evaluate the metric for this point
          vnl_vector<double> xTry = search_fun.GetCoefficients(tSearchTry);
          double f = 0.0;
          search_fun.compute(xTry, &f, NULL);

          // Is this an improvement?
          if(f < fBest)
            {
            fBest = f;
            tLevel->SetMatrix(tSearchTry->GetMatrix());
            tLevel->SetOffset(tSearchTry->GetOffset());
            std::cout << "Rigid search -> Iter " << i << ": " << fBest << " "
                      << xTry << " det = " << vnl_determinant(Qp_search)
                      <<  std::endl;
            }
          }
        }
      }
    else
      {
      // Update the transform from the last level
      MapPhysicalRASSpaceToAffine(of_helper, level, Q_physical, tLevel);
      }

    // Test derivatives
    // Convert to a parameter vector
    vnl_vector<double> xLevel = acf->GetCoefficients(tLevel.GetPointer());

    if(param.flag_debug_deriv)
      {
      // Test the gradient computation
      vnl_vector<double> xGrad(acf->get_number_of_unknowns(), 0.0);
      double f0;
      acf->compute(xLevel, &f0, &xGrad);

      // Propagate the jitter to the transform
      Q_physical = MapAffineToPhysicalRASSpace(of_helper, level, tLevel);
      std::cout << "Initial RAS Transform: " << std::endl << Q_physical  << std::endl;

      printf("ANL gradient: ");
      for (unsigned i = 0; i < xGrad.size(); i++)
        printf("%11.4f ", xGrad[i]);
      printf("\n");

      vnl_vector<double> xGradN(acf->get_number_of_unknowns(), 0.0);
      for(int i = 0; i < acf->get_number_of_unknowns(); i++)
        {
        // double eps = (i % VDim == 0) ? 1.0e-2 : 1.0e-5;
        double eps = param.deriv_epsilon;
        double f1, f2, f3, f4;
        vnl_vector<double> x1 = xLevel, x2 = xLevel, x3 = xLevel, x4 = xLevel;
        x1[i] -= 2 * eps; x2[i] -= eps; x3[i] += eps; x4[i] += 2 * eps;

        // Keep track of gradient even though we do not need it. There is an apparent bug
        // at least with the NCC metric, where the reuse of the working image in a scenario
        // where you first compute the gradient and then do not, causes the iteration through
        // the working image to incorrectly align the per-pixel arrays. Asking for gradient
        // every time is a little more costly, but it avoids this issue
        vnl_vector<double> xGradDummy(acf->get_number_of_unknowns(), 0.0);

        // Four-point derivative computation        
        acf->compute(x1, &f1, &xGradDummy);
        acf->compute(x2, &f2, &xGradDummy);
        acf->compute(x3, &f3, &xGradDummy);
        acf->compute(x4, &f4, &xGradDummy);

        xGradN[i] = (f1 - 8 * f2 + 8 * f3 - f4) / (12 * eps);
        }

      printf("NUM gradient: ");
      for (unsigned i = 0; i < xGradN.size(); i++)
        printf("%11.4f ", xGradN[i]);
      printf("\n");

      std::cout << "f = " << f0 << std::endl;

      acf->GetTransform(xGrad, tLevel.GetPointer());
      std::cout << "A: " << std::endl
                << tLevel->GetMatrix() << std::endl
                << tLevel->GetOffset() << std::endl;

      acf->GetTransform(xGradN, tLevel.GetPointer());
      std::cout << "N: " << std::endl
                << tLevel->GetMatrix() << std::endl
                << tLevel->GetOffset() << std::endl;
      }

    if(param.flag_debug_aff_obj)
      {
      for(int k = -50; k < 50; k++)
        {
        printf("Obj\t%d\t", k);
        for(int i = 0; i < acf->get_number_of_unknowns(); i++)
          {
          vnl_vector<double> xTest = xLevel;
          xTest[i] = xLevel[i] + k * param.deriv_epsilon;
          double f; acf->compute(xTest, &f, NULL);
          printf("%12.8f\t", f);
          }
        printf("\n");
        }
        {
        vnl_vector<double> xTest = xLevel;
          {
          }
        printf("\n");
        }
      }

    // Run the minimization
    if(param.iter_per_level[level] > 0)
      {
      if(param.flag_powell)
        {
        // Set up the optimizer
        vnl_powell *optimizer = new vnl_powell(acf);
        optimizer->set_f_tolerance(1e-9);
        optimizer->set_x_tolerance(1e-4);
        optimizer->set_g_tolerance(1e-6);
        optimizer->set_trace(true);
        optimizer->set_verbose(true);
        optimizer->set_max_function_evals(param.iter_per_level[level]);

        optimizer->minimize(xLevel);
        delete optimizer;

        }
      else
        {
        // Set up the optimizer
        vnl_lbfgs *optimizer = new vnl_lbfgs(*acf);
        optimizer->set_f_tolerance(1e-9);
        optimizer->set_x_tolerance(1e-4);
        optimizer->set_g_tolerance(1e-6);
        optimizer->set_trace(true);
        optimizer->set_max_function_evals(param.iter_per_level[level]);

        optimizer->minimize(xLevel);
        delete optimizer;
        }

      // Did the registration succeed?
      if(xLevel.size() > 0)
        {
        // Get the final transform
        typename LinearTransformType::Pointer tFinal = LinearTransformType::New();
        acf->GetTransform(xLevel, tFinal.GetPointer());
        Q_physical = MapAffineToPhysicalRASSpace(of_helper, level, tFinal);
        }
      else
        {
        // Use the pre-initialization transform parameters
        Q_physical = MapAffineToPhysicalRASSpace(of_helper, level, tLevel);
        }

      // End of level report
      printf("END OF LEVEL %3d\n", level);

      // Print final metric report
      MultiComponentMetricReport metric_report = this->GetMetricLog()[level].back();
      printf("Level %3d  LastIter   Metrics", level);
      for (unsigned i = 0; i < metric_report.ComponentMetrics.size(); i++)
        printf("  %8.6f", metric_report.ComponentMetrics[i]);
      printf("  Energy = %8.6f\n", metric_report.TotalMetric);
      fflush(stdout);
      }

    // Print the final RAS transform for this level (even if no iter)
    printf("Level %3d  Final RAS Transform:\n", level);
    for(unsigned int a = 0; a < VDim+1; a++)
      {
      for(unsigned int b = 0; b < VDim+1; b++)
        printf("%8.4f%c", Q_physical(a,b), b < VDim ? ' ' : '\n');
      }

    delete acf;
    delete pure_acf;
    }

  // Write the final affine transform
  this->WriteAffineMatrixViaCache(param.output, Q_physical);
  return 0;
}





#include "itkStatisticsImageFilter.h"

/** My own time probe because itk's use of fork is messing up my debugging */
class GreedyTimeProbe
{
public:
  GreedyTimeProbe();
  void Start();
  void Stop();
  double GetMean() const;
protected:
  double m_TotalTime;
  double m_StartTime;
  unsigned long m_Runs;
};

GreedyTimeProbe::GreedyTimeProbe()
{
  m_TotalTime = 0.0;
  m_StartTime = 0.0;
  m_Runs = 0.0;
}

void GreedyTimeProbe::Start()
{
  m_StartTime = clock();
}

void GreedyTimeProbe::Stop()
{
  if(m_StartTime == 0.0)
    throw GreedyException("Timer stop without start");
  m_TotalTime += clock() - m_StartTime;
  m_StartTime = 0.0;
  m_Runs++;
}

double GreedyTimeProbe::GetMean() const
{
  if(m_Runs == 0)
    return 0.0;
  else
    return m_TotalTime / (CLOCKS_PER_SEC * m_Runs);
}


template <unsigned int VDim, typename TReal>
std::string
GreedyApproach<VDim, TReal>
::PrintIter(int level, int iter, const MultiComponentMetricReport &metric) const
{
  // Start with a buffer
  char b_level[64], b_iter[64], b_metrics[512], b_line[1024];

  if(level < 0)
    sprintf(b_level, "LastLevel");
  else
    sprintf(b_level, "Level %03d", level);

  if(iter < 0)
    sprintf(b_iter, "LastIter");
  else
    sprintf(b_iter, "Iter %05d", iter);

  if(metric.ComponentMetrics.size() > 1)
    {
    int pos = sprintf(b_metrics, "Metrics");
    for (unsigned i = 0; i < metric.ComponentMetrics.size(); i++)
      pos += sprintf(b_metrics + pos, "  %8.6f", metric.ComponentMetrics[i]);
    }
  else
    sprintf(b_metrics, "");

  sprintf(b_line, "%s  %s  %s  Energy = %8.6f", b_level, b_iter, b_metrics, metric.TotalMetric);
  std::string result = b_line;

  return b_line;
}

/**
 * This is the main function of the GreedyApproach algorithm
 */
template <unsigned int VDim, typename TReal>
int GreedyApproach<VDim, TReal>
::RunDeformable(GreedyParameters &param)
{
  // Create an optical flow helper object
  OFHelperType of_helper;

  // Set the scaling factors for multi-resolution
  of_helper.SetDefaultPyramidFactors(param.iter_per_level.size());

  // Set the scaling mode depending on the metric
  if(param.metric == GreedyParameters::MAHALANOBIS)
    of_helper.SetScaleFixedImageWithVoxelSize(true);

  // Read the image pairs to register
  ReadImages(param, of_helper);

  // An image pointer desribing the current estimate of the deformation
  VectorImagePointer uLevel = NULL;

  // The number of resolution levels
  unsigned nlevels = param.iter_per_level.size();

  // Clear the metric log
  m_MetricLog.clear();

  // Iterate over the resolution levels
  for(unsigned int level = 0; level < nlevels; ++level)
    {
    // Add stage to metric log
    m_MetricLog.push_back(std::vector<MultiComponentMetricReport>());

    // Reference space
    ImageBaseType *refspace = of_helper.GetReferenceSpace(level);

    // Smoothing factors for this level, in physical units
    typename LDDMMType::Vec sigma_pre_phys =
        of_helper.GetSmoothingSigmasInPhysicalUnits(level, param.sigma_pre.sigma,
                                                    param.sigma_pre.physical_units);

    typename LDDMMType::Vec sigma_post_phys =
        of_helper.GetSmoothingSigmasInPhysicalUnits(level, param.sigma_post.sigma,
                                                    param.sigma_post.physical_units);

    // Report the smoothing factors used
    std::cout << "LEVEL " << level+1 << " of " << nlevels << std::endl;
    std::cout << "  Smoothing sigmas: " << sigma_pre_phys << ", " << sigma_post_phys << std::endl;

    // Set up timers for different critical components of the optimization
    GreedyTimeProbe tm_Gradient, tm_Gaussian1, tm_Gaussian2, tm_Iteration,
      tm_Integration, tm_Update;

    // Intermediate images
    ImagePointer iTemp = ImageType::New();
    VectorImagePointer viTemp = VectorImageType::New();
    VectorImagePointer uk = VectorImageType::New();
    VectorImagePointer uk1 = VectorImageType::New();

    // This is the exponentiated uk, in stationary velocity mode it is uk^(2^N)
    VectorImagePointer uk_exp = VectorImageType::New();

    // A pointer to the full warp image - either uk in greedy mode, or uk_exp in diff demons mdoe
    VectorImageType *uFull;

    // Matrix work image (for Lie Bracket) 
    typedef typename LDDMMType::MatrixImageType MatrixImageType;
    typename MatrixImageType::Pointer work_mat = MatrixImageType::New();

    // Allocate the intermediate data
    LDDMMType::alloc_vimg(uk, refspace);
    if(param.iter_per_level[level] > 0)
      {
      LDDMMType::alloc_img(iTemp, refspace);
      LDDMMType::alloc_vimg(viTemp, refspace);
      LDDMMType::alloc_vimg(uk1, refspace);

      // These are only allocated in diffeomorphic demons mode
      if(param.flag_stationary_velocity_mode)
        {
        LDDMMType::alloc_vimg(uk_exp, refspace);
        LDDMMType::alloc_mimg(work_mat, refspace);
        }
      }

    // Initialize the deformation field from last iteration
    if(uLevel.IsNotNull())
      {
      LDDMMType::vimg_resample_identity(uLevel, refspace, uk);
      LDDMMType::vimg_scale_in_place(uk, 2.0);
      uLevel = uk;
      }
    else if(param.initial_warp.size())
      {
      // The user supplied an initial warp or initial root warp. In this case, we
      // do not start iteration from zero, but use the initial warp to start from
      VectorImagePointer uInit = VectorImageType::New();

      // Read the warp file
      LDDMMType::vimg_read(param.initial_warp.c_str(), uInit );

      // Convert the warp file into voxel units from physical units
      OFHelperType::PhysicalWarpToVoxelWarp(uInit, uInit, uInit);

      // Scale the initial warp by the pyramid level
      LDDMMType::vimg_resample_identity(uInit, refspace, uk);
      LDDMMType::vimg_scale_in_place(uk, 1.0 / (1 << level));
      uLevel = uk;
      itk::Index<VDim> test; test.Fill(24);
      std::cout << "Index 24x24x24 maps to " << uInit->GetPixel(test) << std::endl;
      std::cout << "Index 24x24x24 maps to " << uk->GetPixel(test) << std::endl;
      }
    else if(param.affine_init_mode != VOX_IDENTITY)
      {
      typename LinearTransformType::Pointer tran = LinearTransformType::New();

      if(param.affine_init_mode == RAS_FILENAME)
        {
        // Read the initial affine transform from a file
        vnl_matrix<double> Qp = ReadAffineMatrixViaCache(param.affine_init_transform);

        // Map this to voxel space
        MapPhysicalRASSpaceToAffine(of_helper, level, Qp, tran);
        }
      else if(param.affine_init_mode == RAS_IDENTITY)
        {
        // Physical space transform
        vnl_matrix<double> Qp(VDim+1, VDim+1); Qp.set_identity();

        // Map this to voxel space
        MapPhysicalRASSpaceToAffine(of_helper, level, Qp, tran);
        }

      // Create an initial warp
      OFHelperType::AffineToField(tran, uk);
      uLevel = uk;

      itk::Index<VDim> test; test.Fill(24);
      std::cout << "Index 24x24x24 maps to " << uk->GetPixel(test) << std::endl;
      }

    if (uLevel.IsNotNull())
      //LDDMMType::vimg_write(uLevel, "/tmp/ulevel.nii.gz");
      std::cout << "uLevel present but not writing image.\n";

    // Iterate for this level
    for(unsigned int iter = 0; iter < param.iter_per_level[level]; iter++)
      {
      // Start the iteration timer
      tm_Iteration.Start();

      // The epsilon for this level
      double eps= param.epsilon_per_level[level];

      // Integrate the total deformation field for this iteration
      if(param.flag_stationary_velocity_mode)
        {
        tm_Integration.Start();

        // This is the exponentiation of the stationary velocity field
        // Take current warp to 'exponent' power - this is the actual warp
        LDDMMType::vimg_exp(uk, uk_exp, viTemp, param.warp_exponent, 1.0);
        uFull = uk_exp;

        tm_Integration.Stop();
        }
      else
        {
        uFull = uk;
        }

      // Create a metric report that will be returned by all metrics
      MultiComponentMetricReport metric_report;

      // Begin gradient computation
      tm_Gradient.Start();

      // Switch based on the metric
      if(param.metric == GreedyParameters::SSD)
        {
        of_helper.ComputeOpticalFlowField(level, uFull, iTemp, metric_report, uk1, eps);
        metric_report.Scale(1.0 / eps);

        // If there is a mask, multiply the gradient by the mask
        if(param.gradient_mask.size())
          LDDMMType::vimg_multiply_in_place(uk1, of_helper.GetGradientMask(level));
        }

      else if(param.metric == GreedyParameters::MI || param.metric == GreedyParameters::NMI)
        {
        of_helper.ComputeMIFlowField(level, param.metric == GreedyParameters::NMI, uFull, iTemp, metric_report, uk1, eps);

        // If there is a mask, multiply the gradient by the mask
        if(param.gradient_mask.size())
          LDDMMType::vimg_multiply_in_place(uk1, of_helper.GetGradientMask(level));
        }

      else if(param.metric == GreedyParameters::NCC)
        {
        itk::Size<VDim> radius = array_caster<VDim>::to_itkSize(param.metric_radius);

        // Compute the metric - no need to multiply by the mask, this happens already in the NCC metric code
        of_helper.ComputeNCCMetricImage(level, uFull, radius, iTemp, metric_report, uk1, eps);
        metric_report.Scale(1.0 / eps);
        }
      else if(param.metric == GreedyParameters::MAHALANOBIS)
        {
        of_helper.ComputeMahalanobisMetricImage(level, uFull, iTemp, metric_report, uk1);
        }

      // End gradient computation
      tm_Gradient.Stop();

      // Print a report for this iteration
      std::cout << this->PrintIter(level, iter, metric_report) << std::endl;
      fflush(stdout);

      // Record the metric value in the log
      this->RecordMetricValue(metric_report);

      // Dump the gradient image if requested
      if(param.flag_dump_moving && 0 == iter % param.dump_frequency)
        {
        char fname[256];
        sprintf(fname, "dump_gradient_lev%02d_iter%04d.nii.gz", level, iter);
        LDDMMType::vimg_write(uk1, fname);
        }

      // We have now computed the gradient vector field. Next, we smooth it
      tm_Gaussian1.Start();
      LDDMMType::vimg_smooth_withborder(uk1, viTemp, sigma_pre_phys, 1);
      tm_Gaussian1.Stop();

      // After smoothing, compute the maximum vector norm and use it as a normalizing
      // factor for the displacement field
      if(param.time_step_mode == GreedyParameters::SCALE)
        LDDMMType::vimg_normalize_to_fixed_max_length(viTemp, iTemp, eps, false);
      else if (param.time_step_mode == GreedyParameters::SCALEDOWN)
        LDDMMType::vimg_normalize_to_fixed_max_length(viTemp, iTemp, eps, true);

      // Dump the smoothed gradient image if requested
      if(param.flag_dump_moving && 0 == iter % param.dump_frequency)
        {
        char fname[256];
        sprintf(fname, "dump_optflow_lev%02d_iter%04d.nii.gz", level, iter);
        LDDMMType::vimg_write(viTemp, fname);
        }

      // Compute the updated deformation field - in uk1
      tm_Update.Start();
      if(param.flag_stationary_velocity_mode)
        {
        // this is diffeomorphic demons - Vercauteren 2008
        // We now hold the update field in viTemp. This update u should be integrated
        // with the current stationary velocity field such that exp[v'] = exp[v] o exp[u]
        // Vercauteren (2008) suggests using the following expressions
        // v' = v + u (so-so)
        // v' = v + u + [v, u]/2 (this is the Lie bracket)
        
        // Scale the update by 1 / 2^exponent (tiny update, first order approximation)
        LDDMMType::vimg_scale_in_place(viTemp, 1.0 / (2 << param.warp_exponent));

        // Use appropriate update
        if(param.flag_stationary_velocity_mode_use_lie_bracket)
          {
          // Use the Lie Bracket approximation (v + u + [v,u])
          LDDMMType::lie_bracket(uk, viTemp, work_mat, uk1);
          LDDMMType::vimg_scale_in_place(uk1, 0.5); 
          LDDMMType::vimg_add_in_place(uk1, uk);
          LDDMMType::vimg_add_in_place(uk1, viTemp);
          }
        else
          {
          LDDMMType::vimg_copy(uk, uk1);
          LDDMMType::vimg_add_in_place(uk1, viTemp);
          }
        }
      else
        {
        // This is compositive (uk1 = viTemp + uk o viTemp), which is what is done with
        // compositive demons and ANTS
        LDDMMType::interp_vimg(uk, viTemp, 1.0, uk1);
        LDDMMType::vimg_add_in_place(uk1, viTemp);
        }
      tm_Update.Stop();

      // Dump if requested
      if(param.flag_dump_moving && 0 == iter % param.dump_frequency)
        {
        char fname[256];
        sprintf(fname, "dump_uk1_lev%02d_iter%04d.nii.gz", level, iter);
        LDDMMType::vimg_write(uk1, fname);
        }

      // Another layer of smoothing
      tm_Gaussian2.Start();
      LDDMMType::vimg_smooth_withborder(uk1, uk, sigma_post_phys, 1);
      tm_Gaussian2.Stop();

      tm_Iteration.Stop();
      }

    // Store the end result
    uLevel = uk;

    // Compute the jacobian of the deformation field - but only if we iterated at this level
    if(param.iter_per_level[level] > 0)
      {
      LDDMMType::field_jacobian_det(uk, iTemp);
      TReal jac_min, jac_max;
      LDDMMType::img_min_max(iTemp, jac_min, jac_max);
      printf("END OF LEVEL %3d    DetJac Range: %8.4f  to %8.4f \n", level, jac_min, jac_max);

      // Print final metric report
      MultiComponentMetricReport metric_report = this->GetMetricLog()[level].back();
      std::cout << this->PrintIter(level, -1, metric_report) << std::endl;
      fflush(stdout);

      // Print timing information
      printf("  Avg. Gradient Time  : %6.4fs  %5.2f%% \n", tm_Gradient.GetMean(), 
             tm_Gradient.GetMean() * 100.0 / tm_Iteration.GetMean());
      printf("  Avg. Gaussian Time  : %6.4fs  %5.2f%% \n", tm_Gaussian1.GetMean() + tm_Gaussian2.GetMean(),
             (tm_Gaussian1.GetMean() + tm_Gaussian2.GetMean()) * 100.0 / tm_Iteration.GetMean());
      printf("  Avg. Integration Time  : %6.4fs  %5.2f%% \n", tm_Integration.GetMean() + tm_Update.GetMean(),
             (tm_Integration.GetMean() + tm_Update.GetMean()) * 100.0 / tm_Iteration.GetMean());
      printf("  Avg. Total Iteration Time : %6.4fs \n", tm_Iteration.GetMean());
      }
    }

  // The transformation field is in voxel units. To work with ANTS, it must be mapped
  // into physical offset units - just scaled by the spacing?
  
  if(param.flag_stationary_velocity_mode)
    {
    // Take current warp to 'exponent' power - this is the actual warp
    VectorImagePointer uLevelExp = LDDMMType::new_vimg(uLevel);
    VectorImagePointer uLevelWork = LDDMMType::new_vimg(uLevel);
    LDDMMType::vimg_exp(uLevel, uLevelExp, uLevelWork, param.warp_exponent, 1.0);

    // Write the resulting transformation field
    of_helper.WriteCompressedWarpInPhysicalSpace(nlevels - 1, uLevelExp, param.output.c_str(), param.warp_precision);

    if(param.root_warp.size())
      {
      // If asked to write root warp, do so
      of_helper.WriteCompressedWarpInPhysicalSpace(nlevels - 1, uLevel, param.root_warp.c_str(), 0);
      }
    if(param.inverse_warp.size())
      {
      // Compute the inverse (this is probably unnecessary for small warps)
      of_helper.ComputeDeformationFieldInverse(uLevel, uLevelWork, 0);
      of_helper.WriteCompressedWarpInPhysicalSpace(nlevels - 1, uLevelWork, param.inverse_warp.c_str(), param.warp_precision);
      }
    }
  else
    {
    // Write the resulting transformation field
    of_helper.WriteCompressedWarpInPhysicalSpace(nlevels - 1, uLevel, param.output.c_str(), param.warp_precision);

    // If an inverse is requested, compute the inverse using the Chen 2008 fixed method.
    // A modification of this method is that if convergence is slow, we take the square
    // root of the forward transform.
    //
    // TODO: it would be more efficient to check the Lipschitz condition rather than
    // the brute force approach below
    //
    // TODO: the maximum checks should only be done over the region where the warp is
    // not going outside of the image. Right now, they are meaningless and we are doing
    // extra work when computing the inverse.
    if(param.inverse_warp.size())
      {
      // Compute the inverse
      VectorImagePointer uInverse = LDDMMType::new_vimg(uLevel);
      of_helper.ComputeDeformationFieldInverse(uLevel, uInverse, param.warp_exponent);

      // Write the warp using compressed format
      of_helper.WriteCompressedWarpInPhysicalSpace(nlevels - 1, uInverse, param.inverse_warp.c_str(), param.warp_precision);
      }
    }
  return 0;
}




/**
 * Computes the metric without running any optimization. Metric can be computed
 * at different levels by specifying the iterations array
 */
template <unsigned int VDim, typename TReal>
int GreedyApproach<VDim, TReal>
::ComputeMetric(GreedyParameters &param, MultiComponentMetricReport &metric_report)
{
  // Create an optical flow helper object
  OFHelperType of_helper;

  // Set the scaling factors for multi-resolution
  of_helper.SetDefaultPyramidFactors(1);

  // Set the scaling mode depending on the metric
  if(param.metric == GreedyParameters::MAHALANOBIS)
    of_helper.SetScaleFixedImageWithVoxelSize(true);

  // Read the image pairs to register
  ReadImages(param, of_helper);

  // An image pointer desribing the current estimate of the deformation
  VectorImagePointer uLevel = NULL;

  // Reference space
  ImageBaseType *refspace = of_helper.GetReferenceSpace(0);

  // Intermediate images
  ImagePointer iTemp = ImageType::New();
  VectorImagePointer viTemp = VectorImageType::New();
  VectorImagePointer uk = VectorImageType::New();
  VectorImagePointer uk1 = VectorImageType::New();

  // This is the exponentiated uk, in stationary velocity mode it is uk^(2^N)
  VectorImagePointer uk_exp = VectorImageType::New();

  // A pointer to the full warp image - either uk in greedy mode, or uk_exp in diff demons mdoe
  VectorImageType *uFull;

  // Allocate the intermediate data
  LDDMMType::alloc_vimg(uk, refspace);
  LDDMMType::alloc_img(iTemp, refspace);
  LDDMMType::alloc_vimg(viTemp, refspace);
  LDDMMType::alloc_vimg(uk1, refspace);

  // These are only allocated in diffeomorphic demons mode
  if(param.flag_stationary_velocity_mode)
    {
    LDDMMType::alloc_vimg(uk_exp, refspace);
    }

  if(param.initial_warp.size())
    {
    // The user supplied an initial warp or initial root warp. In this case, we
    // do not start iteration from zero, but use the initial warp to start from
    VectorImagePointer uInit = VectorImageType::New();

    // Read the warp file
    LDDMMType::vimg_read(param.initial_warp.c_str(), uInit );

    // Convert the warp file into voxel units from physical units
    OFHelperType::PhysicalWarpToVoxelWarp(uInit, uInit, uInit);

    // Scale the initial warp by the pyramid level
    LDDMMType::vimg_resample_identity(uInit, refspace, uk);
    uLevel = uk;
    }
  else if(param.affine_init_mode != VOX_IDENTITY)
    {
    typename LinearTransformType::Pointer tran = LinearTransformType::New();

    if(param.affine_init_mode == RAS_FILENAME)
      {
      // Read the initial affine transform from a file
      vnl_matrix<double> Qp = ReadAffineMatrixViaCache(param.affine_init_transform);

      // Map this to voxel space
      MapPhysicalRASSpaceToAffine(of_helper, 0, Qp, tran);
      }
    else if(param.affine_init_mode == RAS_IDENTITY)
      {
      // Physical space transform
      vnl_matrix<double> Qp(VDim+1, VDim+1); Qp.set_identity();

      // Map this to voxel space
      MapPhysicalRASSpaceToAffine(of_helper, 0, Qp, tran);
      }

    // Create an initial warp
    OFHelperType::AffineToField(tran, uk);
    uLevel = uk;
    }

  // Integrate the total deformation field for this iteration
  if(param.flag_stationary_velocity_mode)
    {
    // This is the exponentiation of the stationary velocity field
    // Take current warp to 'exponent' power - this is the actual warp
    LDDMMType::vimg_exp(uk, uk_exp, viTemp, param.warp_exponent, 1.0);
    uFull = uk_exp;
    }
  else
    {
    uFull = uk;
    }

  // Switch based on the metric
  if(param.metric == GreedyParameters::SSD)
    {
    of_helper.ComputeOpticalFlowField(0, uFull, iTemp, metric_report, uk1, 1.0);
    }
  else if(param.metric == GreedyParameters::MI || param.metric == GreedyParameters::NMI)
    {
    of_helper.ComputeMIFlowField(0, param.metric == GreedyParameters::NMI, uFull, iTemp, metric_report, uk1, 1.0);

    // If there is a mask, multiply the gradient by the mask
    }

  else if(param.metric == GreedyParameters::NCC)
    {
    itk::Size<VDim> radius = array_caster<VDim>::to_itkSize(param.metric_radius);

    // Compute the metric - no need to multiply by the mask, this happens already in the NCC metric code
    of_helper.ComputeNCCMetricImage(0, uFull, radius, iTemp, metric_report, uk1, 1.0);
    }
  else if(param.metric == GreedyParameters::MAHALANOBIS)
    {
    of_helper.ComputeMahalanobisMetricImage(0, uFull, iTemp, metric_report, uk1);
    }

  return 0;
}


/**
 * This function performs brute force search for similar patches. It generates a discrete displacement
 * field where every pixel in the fixed image is matched to the most similar pixel in the moving image
 * within a certain radius
 */
template <unsigned int VDim, typename TReal>
int GreedyApproach<VDim, TReal>
::RunBrute(GreedyParameters &param)
{
  // Check for valid parameters
  if(param.metric != GreedyParameters::NCC)
    {
    std::cerr << "Brute force search requires NCC metric only" << std::endl;
    return -1;
    }

  if(param.brute_search_radius.size() != VDim)
    {
    std::cerr << "Brute force search radius must be same dimension as the images" << std::endl;
    return -1;
    }

  // Create an optical flow helper object
  OFHelperType of_helper;

  // No multi-resolution
  of_helper.SetDefaultPyramidFactors(1);

  // Read the image pairs to register
  ReadImages(param, of_helper);

  // Reference space
  ImageBaseType *refspace = of_helper.GetReferenceSpace(0);

  // Intermediate images
  VectorImagePointer u_best = LDDMMType::new_vimg(refspace);
  VectorImagePointer u_curr = LDDMMType::new_vimg(refspace);
  ImagePointer m_curr = LDDMMType::new_img(refspace);
  ImagePointer m_best = LDDMMType::new_img(refspace);

  // Allocate m_best to a negative value
  m_best->FillBuffer(-100.0);

  // Create a neighborhood for computing offsets
  itk::Neighborhood<float, VDim> dummy_nbr;
  itk::Size<VDim> search_rad = array_caster<VDim>::to_itkSize(param.brute_search_radius);
  itk::Size<VDim> metric_rad = array_caster<VDim>::to_itkSize(param.metric_radius);
  dummy_nbr.SetRadius(search_rad);

  // Iterate over all offsets
  for(int k = 0; k < dummy_nbr.Size(); k++)
    {
    // Get the offset corresponding to this iteration
    itk::Offset<VDim> offset = dummy_nbr.GetOffset(k);

    // Fill the deformation field with this offset
    typename LDDMMType::Vec vec_offset;
    for(int i = 0; i < VDim; i++)
      vec_offset[i] = offset[i];
    u_curr->FillBuffer(vec_offset);

    // Perform interpolation and metric computation
    MultiComponentMetricReport metric_report;
    of_helper.ComputeNCCMetricImage(0, u_curr, metric_rad, m_curr, metric_report);

    // Temp: keep track of number of updates
    unsigned long n_updates = 0;

    // Out of laziness, just take a quick pass over the images
    typename VectorImageType::RegionType rgn = refspace->GetBufferedRegion();
    itk::ImageRegionIterator<VectorImageType> it_u(u_best, rgn);
    itk::ImageRegionConstIterator<ImageType> it_m_curr(m_curr, rgn);
    itk::ImageRegionIterator<ImageType> it_m_best(m_best, rgn);
    for(; !it_m_best.IsAtEnd(); ++it_m_best, ++it_m_curr, ++it_u)
      {
      float v_curr = it_m_curr.Value();
      if(v_curr > it_m_best.Value())
        {
        it_m_best.Set(v_curr);
        it_u.Set(vec_offset);
        ++n_updates;
        }
      }

    std::cout << "offset: " << offset << "     updates: " << n_updates << std::endl;
    }

  LDDMMType::vimg_write(u_best, param.output.c_str());
  LDDMMType::img_write(m_best, "mbest.nii.gz");

  return 0;
}


#include "itkWarpVectorImageFilter.h"
#include "itkWarpImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"


template <unsigned int VDim, typename TReal>
void GreedyApproach<VDim, TReal>
::ReadTransformChain(const std::vector<TransformSpec> &tran_chain,
                     ImageBaseType *ref_space,
                     VectorImagePointer &out_warp)
{
  // Create the initial transform and set it to zero
  out_warp = VectorImageType::New();
  LDDMMType::alloc_vimg(out_warp, ref_space);

  // Read the sequence of transforms
  for(int i = 0; i < tran_chain.size(); i++)
    {
    // Read the next parameter
    std::string tran = tran_chain[i].filename;

    // Determine if it's an affine transform
    if(itk::ImageIOFactory::CreateImageIO(tran.c_str(), itk::ImageIOFactory::ReadMode))
      {
      // Create a temporary warp
      VectorImagePointer warp_tmp = LDDMMType::new_vimg(ref_space);

      // Read the next warp
      VectorImagePointer warp_i = VectorImageType::New();
      LDDMMType::vimg_read(tran.c_str(), warp_i);

      // If there is an exponent on the transform spec, handle it
      if(tran_chain[i].exponent != 1)
        {
        // The exponent may be specified as a negative number, in which case we take the negative
        // input and exponentiate it
        double absexp = fabs(tran_chain[i].exponent);
        double n_real = log(absexp) / log(2.0);
        int n = (int) (n_real + 0.5);
        if(fabs(n - n_real) > 1.0e-4) 
          throw GreedyException("Currently only power of two exponents are supported for warps");

        // Bring the transform into voxel space
        VectorImagePointer warp_exp = LDDMMType::new_vimg(warp_i);
        OFHelperType::PhysicalWarpToVoxelWarp(warp_i, warp_i, warp_i);

        // Square the transform N times (in its own space)
        LDDMMType::vimg_exp(warp_i, warp_exp, warp_tmp, n, tran_chain[i].exponent / absexp);

        // Bring the transform back into physical space
        OFHelperType::VoxelWarpToPhysicalWarp(warp_exp, warp_i, warp_i);
        }

      // Now we need to compose the current transform and the overall warp.
      LDDMMType::interp_vimg(warp_i, out_warp, 1.0, warp_tmp, false, true);
      LDDMMType::vimg_add_in_place(out_warp, warp_tmp);
      }
    else
      {
      // Read the transform as a matrix
      vnl_matrix<double> mat = ReadAffineMatrixViaCache(tran_chain[i]);
      vnl_matrix<double>  A = mat.extract(VDim, VDim);
      vnl_vector<double> b = mat.get_column(VDim).extract(VDim), q;

      // TODO: stick this in a filter to take advantage of threading!
      typedef itk::ImageRegionIteratorWithIndex<VectorImageType> IterType;
      for(IterType it(out_warp, out_warp->GetBufferedRegion()); !it.IsAtEnd(); ++it)
        {
        itk::Point<double, VDim> pt, pt2;
        typename VectorImageType::IndexType idx = it.GetIndex();

        // Get the physical position
        // TODO: this calls IsInside() internally, which limits efficiency
        out_warp->TransformIndexToPhysicalPoint(idx, pt);

        // Add the displacement (in DICOM coordinates) and
        for(int i = 0; i < VDim; i++)
          pt2[i] = pt[i] + it.Value()[i];

        // Switch to NIFTI coordinates
        pt2[0] = -pt2[0]; pt2[1] = -pt2[1];

        // Apply the matrix - get the transformed coordinate in DICOM space
        q = A * pt2.GetVnlVector() + b;
        q[0] = -q[0]; q[1] = -q[1];

        // Compute the difference in DICOM space
        for(int i = 0; i < VDim; i++)
          it.Value()[i] = q[i] - pt[i];
        }
      }
    }
}

#include "itkBinaryThresholdImageFilter.h"
//#include "itkRecursiveGaussianImageFilter.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkNaryFunctorImageFilter.h"

template <class TInputImage, class TOutputImage>
class NaryLabelVotingFunctor
{
public:
  typedef NaryLabelVotingFunctor<TInputImage,TOutputImage> Self;
  typedef typename TInputImage::PixelType InputPixelType;
  typedef typename TOutputImage::PixelType OutputPixelType;
  typedef std::vector<OutputPixelType> LabelArray;

  NaryLabelVotingFunctor(const LabelArray &labels)
    : m_LabelArray(labels), m_Size(labels.size()) {}

  NaryLabelVotingFunctor() : m_Size(0) {}


  OutputPixelType operator() (const std::vector<InputPixelType> &pix)
  {
    InputPixelType best_val = pix[0];
    int best_index = 0;
    for(int i = 1; i < m_Size; i++)
      if(pix[i] > best_val)
        {
        best_val = pix[i];
        best_index = i;
        }

    return m_LabelArray[best_index];
  }

  bool operator != (const Self &other)
    { return other.m_LabelArray != m_LabelArray; }

protected:
  LabelArray m_LabelArray;
  int m_Size;
};

#include "itkMeshFileReader.h"
#include "itkMeshFileWriter.h"
#include "itkMesh.h"
#include "itkTransformMeshFilter.h"

template <unsigned int VDim, typename TArray>
class PhysicalCoordinateTransform
{
  static void ras_to_lps(const TArray &src, TArray &trg) {}
  static void lps_to_ras(const TArray &src, TArray &trg) {}
};

template <typename TArray>
class PhysicalCoordinateTransform<2, TArray>
{
public:
  static void ras_to_lps(const TArray &src, TArray &trg)
  {
    trg[0] = -src[0];
    trg[1] = -src[1];
  }

  static void lps_to_ras(const TArray &src, TArray &trg)
  {
    trg[0] = -src[0];
    trg[1] = -src[1];
  }
};

template <typename TArray>
class PhysicalCoordinateTransform<3, TArray>
{
public:
  static void ras_to_lps(const TArray &src, TArray &trg)
  {
    trg[0] = -src[0];
    trg[1] = -src[1];
    trg[2] = src[2];
  }

  static void lps_to_ras(const TArray &src, TArray &trg)
  {
    trg[0] = -src[0];
    trg[1] = -src[1];
    trg[2] = src[2];
  }
};

template <typename TArray>
class PhysicalCoordinateTransform<4, TArray>
{
public:
  static void ras_to_lps(const TArray &src, TArray &trg)
  {
    trg[0] = -src[0];
    trg[1] = -src[1];
    trg[2] = src[2];
    trg[3] = src[3];
  }

  static void lps_to_ras(const TArray &src, TArray &trg)
  {
    trg[0] = -src[0];
    trg[1] = -src[1];
    trg[2] = src[2];
    trg[3] = src[3];
  }
};



template <unsigned int VDim, typename TReal>
class WarpMeshTransformFunctor : public itk::DataObject
{
public:
  typedef WarpMeshTransformFunctor<VDim, TReal>       Self;
  typedef itk::DataObject                             Superclass;
  typedef itk::SmartPointer<Self>                     Pointer;
  typedef itk::SmartPointer<const Self>               ConstPointer;

  itkTypeMacro(WarpMeshTransformFunctor, itk::DataObject)
  itkNewMacro(Self)

  typedef GreedyApproach<VDim, TReal> GreedyAPI;
  typedef typename GreedyAPI::VectorImageType VectorImageType;
  typedef typename GreedyAPI::ImageBaseType ImageBaseType;
  typedef FastLinearInterpolator<VectorImageType, TReal, VDim> FastInterpolator;
  typedef itk::ContinuousIndex<TReal, VDim> CIndexType;
  typedef itk::Point<TReal, VDim> PointType;

  void SetWarp(VectorImageType *warp)
  {
    if(m_Interpolator) delete m_Interpolator;
    m_Interpolator = new FastInterpolator(warp);
  }

  void SetReferenceSpace(ImageBaseType *ref)
  {
    m_ReferenceSpace = ref;
  }

  PointType TransformPoint(const PointType &x)
  {
    // Our convention is to use NIFTI/RAS coordinates for meshes, whereas ITK
    // uses the DICOM/LPS convention. We transform point to LPS first
    PointType x_lps, phi_x;

    PhysicalCoordinateTransform<VDim, PointType>::ras_to_lps(x, x_lps);

    CIndexType cix;
    typename VectorImageType::PixelType vec;
    vec.Fill(0.0);
    m_ReferenceSpace->TransformPhysicalPointToContinuousIndex(x_lps, cix);
    m_Interpolator->Interpolate(cix.GetDataPointer(), &vec);

    for(int d = 0; d < VDim; d++)
      {
      phi_x[d] = vec[d] + x_lps[d];
      }


    PhysicalCoordinateTransform<VDim, PointType>::lps_to_ras(phi_x, phi_x);

    return phi_x;
  }

protected:

  WarpMeshTransformFunctor() { m_Interpolator = NULL; }
  ~WarpMeshTransformFunctor()
  {
    if(m_Interpolator)
      delete m_Interpolator;
  }

private:

  typename ImageBaseType::Pointer m_ReferenceSpace;
  FastInterpolator *m_Interpolator;

};

/**
 * This code computes the jacobian determinant field for a deformation. The
 * recommended mode for this computation is to take the k-th root of the input
 * transformation and then compose the Jacobians
 */
template <unsigned int VDim, typename TReal>
int GreedyApproach<VDim, TReal>
::RunJacobian(GreedyParameters &param)
{
  // Read the warp as a transform chain
  VectorImagePointer warp;

  // Read the warp file
  LDDMMType::vimg_read(param.jacobian_param.in_warp.c_str(), warp);

  // Convert the warp file into voxel units from physical units
  OFHelperType::PhysicalWarpToVoxelWarp(warp, warp, warp);

  // Compute the root of the warp
  VectorImagePointer root_warp = VectorImageType::New();
  LDDMMType::alloc_vimg(root_warp, warp);

  // Allocate a working warp
  VectorImagePointer work_warp = VectorImageType::New();
  LDDMMType::alloc_vimg(work_warp, warp);

  // Compute the root warp, which is not stored in the variable warp
  OFHelperType::ComputeWarpRoot(warp, root_warp, param.warp_exponent);

  // Initialize empty array of Jacobians
  typedef typename LDDMMType::MatrixImageType JacobianImageType;
  typename JacobianImageType::Pointer jac = LDDMMType::new_mimg(warp);

  typename JacobianImageType::Pointer jac_work = LDDMMType::new_mimg(warp);

  // Compute the Jacobian of the root warp
  LDDMMType::field_jacobian(root_warp, jac);

  // Compute the Jacobian matrix of the root warp; jac[a] = D_a (warp)
  for(int k = 0; k < param.warp_exponent; k++)
    {
    // Compute the composition of the Jacobian with itself
    LDDMMType::jacobian_of_composition(jac, jac, root_warp, jac_work);

    // Swap the pointers, so jac points to the actual composed jacobian
    typename JacobianImageType::Pointer temp = jac_work.GetPointer();
    jac_work = jac.GetPointer();
    jac = temp.GetPointer();

    // Compute the composition of the warp with itself, place into root_warp
    LDDMMType::interp_vimg(root_warp, root_warp, 1.0, work_warp);
    LDDMMType::vimg_add_in_place(root_warp, work_warp);
    }

  // At this point, root_warp should hold the original warp, and jac+I will hold
  // the Jacobian of the original warp. We need to compute the determinant
  ImagePointer jac_det = ImageType::New();
  LDDMMType::alloc_img(jac_det, warp);
  LDDMMType::mimg_det(jac, 1.0, jac_det);

  // Write the computed Jacobian
  LDDMMType::img_write(jac_det, param.jacobian_param.out_det_jac.c_str(), itk::ImageIOBase::FLOAT);
  return 0;
}

/**
 * Run the reslice code - simply apply a warp or set of warps to images
 */
template <unsigned int VDim, typename TReal>
int GreedyApproach<VDim, TReal>
::RunReslice(GreedyParameters &param)
{
  GreedyResliceParameters r_param = param.reslice_param;

  // Check the parameters
  if(!r_param.ref_image.size())
    throw GreedyException("A reference image (-rf) option is required for reslice commands");

  if(r_param.images.size() + r_param.meshes.size() == 0
     && !r_param.out_composed_warp.size()
     && !r_param.out_jacobian_image.size())
    throw GreedyException("No operation specified for reslice mode. "
                          "Use one of -rm, -rs or -rc commands.");

  // Read the fixed as a plain image (we don't care if it's composite)
  typename ImageBaseType::Pointer ref = ReadImageBaseViaCache(r_param.ref_image);

  // Read the transform chain
  VectorImagePointer warp;
  ReadTransformChain(param.reslice_param.transforms, ref, warp);

  // Write the composite warp if requested
  if(r_param.out_composed_warp.size())
    {
    WriteImageViaCache(warp.GetPointer(), r_param.out_composed_warp.c_str(), itk::ImageIOBase::FLOAT);
    }

  // Compute the Jacobian of the warp if requested
  if(r_param.out_jacobian_image.size())
    {
    ImagePointer iTemp = ImageType::New();
    LDDMMType::alloc_img(iTemp, warp);
    LDDMMType::field_jacobian_det(warp, iTemp);

    WriteImageViaCache(iTemp.GetPointer(), r_param.out_jacobian_image.c_str(), itk::ImageIOBase::FLOAT);
    }


  // Process image pairs
  for(int i = 0; i < r_param.images.size(); i++)
    {
    const char *filename = r_param.images[i].moving.c_str();

    // Handle the special case of multi-label images
    if(r_param.images[i].interp.mode == InterpSpec::LABELWISE)
      {
      // The label image assumed to be an image of shortsC
      typedef itk::Image<short, VDim> LabelImageType;
      typedef itk::ImageFileReader<LabelImageType> LabelReaderType;

      // Create a reader
      typename LabelReaderType::Pointer reader = LabelReaderType::New();
      reader->SetFileName(filename);
      reader->Update();
      typename LabelImageType::Pointer moving = reader->GetOutput();

      // Scan the unique labels in the image
      std::set<short> label_set;
      short *labels = moving->GetBufferPointer();
      int n_pixels = moving->GetPixelContainer()->Size();

      // Get the list of unique pixels
      short last_pixel = 0;
      for(int j = 0; j < n_pixels; j++)
        {
        short pixel = labels[j];
        if(last_pixel != pixel || i == 0)
          {
          label_set.insert(pixel);
          last_pixel = pixel;
          if(label_set.size() > 1000)
            throw GreedyException("Label wise interpolation not supported for image %s "
                                  "which has over 1000 distinct labels", filename);
          }
        }

      // Turn this set into an array
      std::vector<short> label_array(label_set.begin(), label_set.end());

      // Create a N-way voting filter
      typedef NaryLabelVotingFunctor<ImageType, LabelImageType> VotingFunctor;
      VotingFunctor vf(label_array);

      typedef itk::NaryFunctorImageFilter<ImageType, LabelImageType, VotingFunctor> VotingFilter;
      typename VotingFilter::Pointer fltVoting = VotingFilter::New();
      fltVoting->SetFunctor(vf);

      // Create a mini-pipeline of streaming filters
      for(int j = 0; j < label_array.size(); j++)
        {
        // Set up a threshold filter for this label
        typedef itk::BinaryThresholdImageFilter<LabelImageType, ImageType> ThresholdFilterType;
        typename ThresholdFilterType::Pointer fltThreshold = ThresholdFilterType::New();
        fltThreshold->SetInput(moving);
        fltThreshold->SetLowerThreshold(label_array[j]);
        fltThreshold->SetUpperThreshold(label_array[j]);
        fltThreshold->SetInsideValue(1.0);
        fltThreshold->SetOutsideValue(0.0);

        // Set up a smoothing filter for this label
        typedef itk::SmoothingRecursiveGaussianImageFilter<ImageType, ImageType> SmootherType;
        typename SmootherType::Pointer fltSmooth = SmootherType::New();
        fltSmooth->SetInput(fltThreshold->GetOutput());

        // Work out the sigmas for the filter
        if(r_param.images[i].interp.sigma.physical_units)
          {
          fltSmooth->SetSigma(r_param.images[i].interp.sigma.sigma);
          }
        else
          {
          typename SmootherType::SigmaArrayType sigma_array;
          for(int d = 0; d < VDim; d++)
            sigma_array[d] = r_param.images[i].interp.sigma.sigma * moving->GetSpacing()[d];
          fltSmooth->SetSigmaArray(sigma_array);
          }

        // TODO: we should really be coercing the output into a vector image to speed up interpolation!
        typedef FastWarpCompositeImageFilter<ImageType, ImageType, VectorImageType> InterpFilter;
        typename InterpFilter::Pointer fltInterp = InterpFilter::New();
        fltInterp->SetMovingImage(fltSmooth->GetOutput());
        fltInterp->SetDeformationField(warp);
        fltInterp->SetUsePhysicalSpace(true);

        fltInterp->Update();

        // Add to the voting filter
        fltVoting->SetInput(j, fltInterp->GetOutput());
        }

      // TODO: test out streaming!
      // Run this big pipeline
      fltVoting->Update();

      // Save
      WriteImageViaCache(fltVoting->GetOutput(), r_param.images[i].output.c_str());
      }
    else
      {
      // Read the input image and record its type
      itk::ImageIOBase::IOComponentType comp;
      CompositeImagePointer moving = ReadImageViaCache<CompositeImageType>(filename, &comp);

      // Allocate the warped image
      CompositeImagePointer warped = LDDMMType::new_cimg(ref, moving->GetNumberOfComponentsPerPixel());

      // Perform the warp
      LDDMMType::interp_cimg(moving, warp, warped,
                             r_param.images[i].interp.mode == InterpSpec::NEAREST,
                             true, r_param.images[i].interp.outside_value);

      // Write, casting to the input component type
      WriteImageViaCache(warped.GetPointer(), r_param.images[i].output.c_str(), comp);
      }
    }

  // Process meshes
  for(int i = 0; i < r_param.meshes.size(); i++)
    {
    typedef itk::Mesh<TReal, VDim> MeshType;
    typedef itk::MeshFileReader<MeshType> MeshReader;
    typename MeshType::Pointer mesh;

    if(itksys::SystemTools::GetFilenameExtension(r_param.meshes[i].fixed) == ".csv")
      {
      mesh = MeshType::New();

      std::ifstream fin(r_param.meshes[i].fixed.c_str());
      std::string f_line, f_token;
      unsigned int n_pts = 0;
      while(std::getline(fin, f_line))
        {
        std::istringstream iss(f_line);
        itk::Point<TReal, VDim> pt;
        for(unsigned int a = 0; a < VDim; a++)
          {
          if(!std::getline(iss, f_token, ','))
            throw GreedyException("Error reading CSV file, line %s", f_line.c_str());
          pt[a] = atof(f_token.c_str());
          }
        mesh->SetPoint(n_pts++, pt);
        }
      }
    else
      {
      typename MeshReader::Pointer reader = MeshReader::New();
      reader->SetFileName(r_param.meshes[i].fixed.c_str());
      reader->Update();
      mesh = reader->GetOutput();
      }

    typedef WarpMeshTransformFunctor<VDim, TReal> TransformType;
    typename TransformType::Pointer transform = TransformType::New();
    transform->SetWarp(warp);
    transform->SetReferenceSpace(ref);

    typedef itk::TransformMeshFilter<MeshType, MeshType, TransformType> FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    filter->SetTransform(transform);
    filter->SetInput(mesh);
    filter->Update();

    mesh = filter->GetOutput();

    if(itksys::SystemTools::GetFilenameExtension(r_param.meshes[i].output) == ".csv")
      {
      std::ofstream out(r_param.meshes[i].output.c_str());
      for(unsigned int i = 0; i < mesh->GetNumberOfPoints(); i++)
        {
        itk::Point<TReal, VDim> pt = mesh->GetPoint(i);
        for(unsigned int a = 0; a < VDim; a++)
          out << pt[a] << (a < VDim-1 ? "," : "\n");
        }
      }
    else
      {
      typedef itk::MeshFileWriter<MeshType> MeshWriter;
      typename MeshWriter::Pointer writer = MeshWriter::New();
      writer->SetInput(mesh);
      writer->SetFileName(r_param.meshes[i].output.c_str());
      writer->Update();
      }
    }



  return 0;
}

template <unsigned int VDim, typename TReal>
void GreedyApproach<VDim, TReal>
::ComputeImageMoments(CompositeImageType *image,
                      const std::vector<double> &weights,
                      VecFx &m1, MatFx &m2)
{
  int n = image->GetNumberOfComponentsPerPixel();
  TReal sum_I = 0.0;
  m1.fill(0.0); m2.fill(0.0);

  typedef itk::ImageRegionConstIteratorWithIndex<CompositeImageType> Iterator;
  for(Iterator it(image, image->GetBufferedRegion()); !it.IsAtEnd(); ++it)
    {
    typedef itk::Point<TReal, VDim> PointType;
    PointType p_lps, p_ras;
    image->TransformIndexToPhysicalPoint(it.GetIndex(), p_lps);
    PhysicalCoordinateTransform<VDim, PointType>::lps_to_ras(p_lps, p_ras);
    VecFx X(p_ras.GetDataPointer());
    MatFx X2 = outer_product(X, X);

    typename CompositeImageType::PixelType pix = it.Get();

    // Just weight the components of intensity by weight vector - this sort of makes sense?
    TReal val = 0.0;
    for(int k = 0; k < n; k++)
      val += weights[k] * pix[k];

    sum_I += val;
    m1 += X * val;
    m2 += X2 * val;
    }

  // Compute the mean and covariance from the sum of squares
  m1 = m1 / sum_I;
  m2 = (m2 - sum_I *  outer_product(m1, m1)) / sum_I;
}

template <unsigned int VDim, typename TReal>
int GreedyApproach<VDim, TReal>
::RunAlignMoments(GreedyParameters &param)
{
  typedef PhysicalSpaceAffineCostFunction<VDim, TReal> PhysicalSpaceAffineCostFunction;

  // Create an optical flow helper object
  OFHelperType of_helper;

  // No multi-resolution
  of_helper.SetDefaultPyramidFactors(1);

  // Read the image pairs to register
  ReadImages(param, of_helper);

  // Compute the moments of intertia for the fixed and moving images. For now
  // this is done in an iterator loop, out of laziness. Should be converted to
  // a filter if this whole moments business proves useful
  VecFx m1f, m1m;
  MatFx m2f, m2m;


  std::cout << "--- MATCHING BY MOMENTS OF ORDER " << param.moments_order << " ---" << std::endl;

  ComputeImageMoments(of_helper.GetFixedComposite(0), of_helper.GetWeights(), m1f, m2f);

  std::cout << "Fixed Mean        : " << m1f << std::endl;
  std::cout << "Fixed Covariance  : " << std::endl << m2f << std::endl;

  ComputeImageMoments(of_helper.GetMovingComposite(0), of_helper.GetWeights(), m1m, m2m);

  std::cout << "Moving Mean       : " << m1m << std::endl;
  std::cout << "Moving Covariance : " << std::endl << m2m << std::endl;

  // This flag forces no rotation, only flip
  if(param.moments_order == 1 || param.flag_moments_id_covariance)
    {
    m2f.set_identity();
    m2m.set_identity();
    }

  // Decompose covariance matrices into eigenvectors and eigenvalues
  vnl_vector<TReal> Df, Dm;
  vnl_matrix<TReal> Vf, Vm;
  vnl_symmetric_eigensystem_compute<TReal>(m2f, Vf, Df);
  vnl_symmetric_eigensystem_compute<TReal>(m2m, Vm, Dm);

  // Create a rigid registration problem
  PhysicalSpaceAffineCostFunction cost_fn(&param, this, 0, &of_helper);

  // The best set of coefficients and the associated match value
  vnl_vector<double> xBest;
  TReal xBestMatch = vnl_numeric_traits<TReal>::maxval;

  // Generate all possible flip matrices
  int n_flip = 1 << VDim;
  for(int k_flip = 0; k_flip < n_flip; k_flip++)
    {
    // If using first moments, ignore all flips, only allow identity
    if(param.moments_order == 1 && k_flip != n_flip - 1)
      continue;

    // Generate the flip matrix
    MatFx F(0.0);
    for(int d = 0; d < VDim; d++)
      F(d,d) = (k_flip & (1 << d)) ? 1 : -1;;

    // Compute the rotation matrix - takes fixed coordinates into moving space
    MatFx R = Vm * F * Vf.transpose();
    VecFx b = m1m - R * m1f;

    vnl_matrix<TReal> A(VDim+1, VDim+1, 0.0);
    A.set_identity();
    A.update(R, 0, 0);
    for(int d= 0 ;d< VDim;d++)
      A(d,VDim) = b[d];

    // Ignore flips with the wrong determinant
    double det_R = vnl_determinant(R);
    if((param.moments_order == 2 && param.moments_flip_determinant == 1 && det_R < 0) ||
       (param.moments_order == 2 && param.moments_flip_determinant == -1 && det_R > 0))
      {
      continue;
      }

    // Generate affine coefficients from the rotation and shift
    vnl_vector<double> x(cost_fn.get_number_of_unknowns());
    flatten_affine_transform(R, b, x.data_block());

    // Compute similarity
    double f = 0.0;
    cost_fn.compute(x, &f, NULL);

    std::cout << "Metric for flip " << F.get_diagonal() << " : " << f << std::endl;

    // Compare
    if(xBestMatch > f || xBest.size() == 0)
      {
      xBestMatch = f;
      xBest = x;
      }
    }

  // Save the best transform
  typename LinearTransformType::Pointer tran = LinearTransformType::New();
  cost_fn.GetTransform(xBest, tran);
  vnl_matrix<double> Q_physical = MapAffineToPhysicalRASSpace(of_helper, 0, tran);
  this->WriteAffineMatrixViaCache(param.output, Q_physical);

  return 0;
}

/**
 * Post-hoc warp inversion - the Achilles heel of non-symmetric registration :(
 */
template <unsigned int VDim, typename TReal>
int GreedyApproach<VDim, TReal>
::RunInvertWarp(GreedyParameters &param)
{
  // Read the warp as a transform chain
  VectorImagePointer warp;

  // Read the warp file
  LDDMMType::vimg_read(param.invwarp_param.in_warp.c_str(), warp);

  // Convert the warp file into voxel units from physical units
  OFHelperType::PhysicalWarpToVoxelWarp(warp, warp, warp);


  // Compute the inverse of the warp
  VectorImagePointer uInverse = VectorImageType::New();
  LDDMMType::alloc_vimg(uInverse, warp);
  OFHelperType::ComputeDeformationFieldInverse(warp, uInverse, param.warp_exponent, true);

  // Write the warp using compressed format
  OFHelperType::WriteCompressedWarpInPhysicalSpace(uInverse, warp, param.invwarp_param.out_warp.c_str(), param.warp_precision);

  return 0;
}

/**
 * Post-hoc warp root
 */
template <unsigned int VDim, typename TReal>
int GreedyApproach<VDim, TReal>
::RunRootWarp(GreedyParameters &param)
{
  // Read the warp as a transform chain
  VectorImagePointer warp;

  // Read the warp file
  LDDMMType::vimg_read(param.warproot_param.in_warp.c_str(), warp);

  // Convert the warp file into voxel units from physical units
  OFHelperType::PhysicalWarpToVoxelWarp(warp, warp, warp);

  // Allocate the root
  VectorImagePointer warp_root = VectorImageType::New();
  LDDMMType::alloc_vimg(warp_root, warp);

  // Take the n-th root
  OFHelperType::ComputeWarpRoot(warp, warp_root, param.warp_exponent, 1e-6);

  // Write the warp using compressed format
  OFHelperType::WriteCompressedWarpInPhysicalSpace(warp_root, warp, param.warproot_param.out_warp.c_str(), param.warp_precision);

  return 0;
}


template <unsigned int VDim, typename TReal>
void GreedyApproach<VDim, TReal>
::AddCachedInputObject(std::string key, itk::Object *object)
{
  m_ImageCache[key].target = object;
  m_ImageCache[key].force_write = false;
}

template <unsigned int VDim, typename TReal>
void GreedyApproach<VDim, TReal>
::AddCachedOutputObject(std::string key, itk::Object *object, bool force_write)
{
  m_ImageCache[key].target = object;
  m_ImageCache[key].force_write = force_write;
}

template <unsigned int VDim, typename TReal>
const typename GreedyApproach<VDim,TReal>::MetricLogType &
GreedyApproach<VDim,TReal>
::GetMetricLog() const
{
  return m_MetricLog;
}

template<unsigned int VDim, typename TReal>
MultiComponentMetricReport GreedyApproach<VDim, TReal>
::GetLastMetricReport() const
{
  // Find last non-empty result
  for(int k = m_MetricLog.size()-1; k >= 0; --k)
    {
    if(m_MetricLog[k].size())
      return m_MetricLog[k].back();
    }

  // If empty, throw exception
  throw GreedyException("Metric log is empty in GetLastMetricValue()");
  return MultiComponentMetricReport();
}

template <unsigned int VDim, typename TReal>
void GreedyApproach<VDim, TReal>
::ConfigThreads(const GreedyParameters &param)
{
  if(param.threads > 0)
    {
    std::cout << "Limiting the number of threads to " << param.threads << std::endl;
    itk::MultiThreader::SetGlobalMaximumNumberOfThreads(param.threads);
    }
  else
    {
    std::cout << "Executing with the default number of threads: " << itk::MultiThreader::GetGlobalDefaultNumberOfThreads() << std::endl;

    }
}

template<unsigned int VDim, typename TReal>
int GreedyApproach<VDim, TReal>
::Run(GreedyParameters &param)
{
  ConfigThreads(param);

  switch(param.mode)
    {
    case GreedyParameters::GREEDY:
      return Self::RunDeformable(param);
    case GreedyParameters::AFFINE:
      return Self::RunAffine(param);
    case GreedyParameters::BRUTE:
      return Self::RunBrute(param);
    case GreedyParameters::MOMENTS:
      return Self::RunAlignMoments(param);
    case GreedyParameters::RESLICE:
      return Self::RunReslice(param);
    case GreedyParameters::INVERT_WARP:
      return Self::RunInvertWarp(param);
    case GreedyParameters::JACOBIAN_WARP:
      return Self::RunJacobian(param);
    case GreedyParameters::ROOT_WARP:
      return Self::RunRootWarp(param);
    }

  return -1;
}




template class GreedyApproach<2, float>;
template class GreedyApproach<3, float>;
template class GreedyApproach<4, float>;
template class GreedyApproach<2, double>;
template class GreedyApproach<3, double>;
template class GreedyApproach<4, double>;
