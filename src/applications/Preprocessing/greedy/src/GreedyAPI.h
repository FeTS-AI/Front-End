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
#ifndef GREEDYAPI_H
#define GREEDYAPI_H

#include "GreedyParameters.h"
#include "GreedyException.h"
#include "MultiComponentMetricReport.h"
#include "lddmm_data.h"
#include "AffineCostFunctions.h"
#include <vnl/vnl_random.h>
#include <map>
#include "itkCommand.h"

template <typename T, unsigned int V> class MultiImageOpticalFlowHelper;

namespace itk {
  template <typename T, unsigned int D1, unsigned int D2> class MatrixOffsetTransformBase;

}

/**
 * This is the top level class for the greedy software. It contains methods
 * for deformable and affine registration.
 */
template <unsigned int VDim, typename TReal = double>
class GreedyApproach
{
public:

  typedef GreedyApproach<VDim, TReal> Self;

  typedef LDDMMData<TReal, VDim> LDDMMType;
  typedef typename LDDMMType::ImageBaseType ImageBaseType;
  typedef typename LDDMMType::ImageType ImageType;
  typedef typename LDDMMType::ImagePointer ImagePointer;
  typedef typename LDDMMType::VectorImageType VectorImageType;
  typedef typename LDDMMType::VectorImagePointer VectorImagePointer;
  typedef typename LDDMMType::CompositeImageType CompositeImageType;
  typedef typename LDDMMType::CompositeImagePointer CompositeImagePointer;

  typedef vnl_vector_fixed<TReal, VDim> VecFx;
  typedef vnl_matrix_fixed<TReal, VDim, VDim> MatFx;

  typedef std::vector< std::vector<MultiComponentMetricReport> > MetricLogType;

  typedef MultiImageOpticalFlowHelper<TReal, VDim> OFHelperType;

  typedef itk::MatrixOffsetTransformBase<TReal, VDim, VDim> LinearTransformType;

  struct ImagePair {
    ImagePointer fixed, moving;
    VectorImagePointer grad_moving;
    double weight;
  };

  static void ConfigThreads(const GreedyParameters &param);

  int Run(GreedyParameters &param);

  int RunDeformable(GreedyParameters &param);

  int RunAffine(GreedyParameters &param);

  int RunBrute(GreedyParameters &param);

  int RunReslice(GreedyParameters &param);

  int RunInvertWarp(GreedyParameters &param);

  int RunRootWarp(GreedyParameters &param);

  int RunAlignMoments(GreedyParameters &param);

  int RunJacobian(GreedyParameters &param);

  int ComputeMetric(GreedyParameters &param, MultiComponentMetricReport &metric_report);

  /**
   * Add an image that is already in memory to the internal cache, and
   * associate it with a filename. This provides a way for images already
   * loaded in memory to be passed in to the Greedy API while using the
   * standard parameter structures.
   *
   * Normally, images such as the fixed image are passed as part of the
   * GreedyParameters object as filenames. For example, we might set
   *
   *   param.inputs[0].fixed = "/tmp/goo.nii.gz";
   *
   * However, if we are linking to the greedy API from another program and
   * already have the fixed image in memory, we can use the cache mechanism
   * instead.
   *
   *   greedyapi.AddCachedInputObject("FIXED-0", myimage);
   *   param.inputs[0].fixed = "FIXED-0";
   *
   * The API will check the cache before loading the image. The type of the
   * object in the cache must match the type of the object expected internally,
   * which is VectorImage for most images. If not, an exception will be
   * thrown.
   *
   * Note that the cache does not use smart pointers to refer to the objects
   * so it's the caller's responsibility to keep the object pointed to while
   * the API is being used.
   */
  void AddCachedInputObject(std::string key, itk::Object *object);

  /**
   * Add an image/matrix to the output cache. This has the same behavior as
   * the input cache, but there is an additional flag as to whether you want
   * to save the output object to the specified filename in addition to writing
   * it to the cached image/matrix. This allows you to both store the result in
   * the cache and write it to a filename specified in the key
   */
  void AddCachedOutputObject(std::string key, itk::Object *object, bool force_write = false);

  /**
   * Get the metric log - values of metric per level. Can be called from
   * callback functions and observers
   */
  const MetricLogType &GetMetricLog() const;

  /** Get the last value of the metric recorded */
  MultiComponentMetricReport GetLastMetricReport() const;

  vnl_matrix<double> ReadAffineMatrixViaCache(const TransformSpec &ts);

  void WriteAffineMatrixViaCache(const std::string &filename, const vnl_matrix<double> &Qp);

  static vnl_matrix<double> ReadAffineMatrix(const TransformSpec &ts);

  static void WriteAffineMatrix(const std::string &filename, const vnl_matrix<double> &Qp);

  static vnl_matrix<double> MapAffineToPhysicalRASSpace(
      OFHelperType &of_helper, int level,
      LinearTransformType *tran);

  static void MapPhysicalRASSpaceToAffine(
      OFHelperType &of_helper, int level,
      vnl_matrix<double> &Qp,
      LinearTransformType *tran);

  void RecordMetricValue(const MultiComponentMetricReport &metric);

  // Helper method to print iteration reports
  std::string PrintIter(int level, int iter, const MultiComponentMetricReport &metric) const;

protected:

  struct CacheEntry {
    itk::Object *target;
    bool force_write;
  };

  typedef std::map<std::string, CacheEntry> ImageCache;
  ImageCache m_ImageCache;

  // A log of metric values used during registration - so metric can be looked up
  // in the callbacks to RunAffine, etc.
  MetricLogType m_MetricLog;

  // This function reads the image from disk, or from a memory location mapped to a
  // string. The first approach is used by the command-line interface, and the second
  // approach is used by the API, allowing images to be passed from other software.
  // An optional second argument is used to store the component type, but only if
  // the image is actually loaded from disk. For cached images, the component type
  // will be unknown.
  template <class TImage>
  itk::SmartPointer<TImage> ReadImageViaCache(const std::string &filename,
                                              itk::ImageIOBase::IOComponentType *comp_type = NULL);

  // This function reads an image base object via cache. It is more permissive than using
  // ReadImageViaCache.
  typename ImageBaseType::Pointer ReadImageBaseViaCache(const std::string &filename);


  // Write an image using the cache
  template <class TImage>
  void WriteImageViaCache(TImage *img, const std::string &filename,
                          typename LDDMMType::IOComponentType comp = itk::ImageIOBase::UNKNOWNCOMPONENTTYPE);

  void ReadImages(GreedyParameters &param, OFHelperType &ofhelper);

  void ReadTransformChain(const std::vector<TransformSpec> &tran_chain,
                          ImageBaseType *ref_space,
                          VectorImagePointer &out_warp);

  // Compute the moments of a composite image (mean and covariance matrix of coordinate weighted by intensity)
  void ComputeImageMoments(CompositeImageType *image, const std::vector<double> &weights, VecFx &m1, MatFx &m2);



  // friend class PureAffineCostFunction<VDim, TReal>;

};

// Little helper functions
template <unsigned int VDim> class array_caster
{
public:
  template <class T> static itk::Size<VDim> to_itkSize(const T &t)
  {
    itk::Size<VDim> sz;
    for(int i = 0; i < VDim; i++)
      sz[i] = t[i];
    return sz;
  }
};


#endif // GREEDYAPI_H
