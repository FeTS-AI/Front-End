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
#ifndef AFFINECOSTFUNCTIONS_H
#define AFFINECOSTFUNCTIONS_H

#include "GreedyAPI.h"
#include <vnl/vnl_cost_function.h>
#include <vnl/vnl_random.h>
#include <vnl/vnl_trace.h>

namespace itk {
  template <typename T, unsigned int D1, unsigned int D2> class MatrixOffsetTransformBase;
}

template <unsigned int VDim, typename TReal> class GreedyApproach;
template <typename T, unsigned int V> class MultiImageOpticalFlowHelper;

/**
 * Parent of all affine/rigid cost functions
 */
template <unsigned int VDim, typename TReal = double>
class AbstractAffineCostFunction : public vnl_cost_function
{
public:
  typedef itk::MatrixOffsetTransformBase<TReal, VDim, VDim> LinearTransformType;
  typedef MultiImageOpticalFlowHelper<TReal, VDim> OFHelperType;
  typedef GreedyApproach<VDim, TReal> ParentType;

  AbstractAffineCostFunction(int n_unknowns) : vnl_cost_function(n_unknowns) {}
  virtual vnl_vector<double> GetCoefficients(LinearTransformType *tran) = 0;
  virtual void GetTransform(const vnl_vector<double> &coeff, LinearTransformType *tran) = 0;
  virtual void compute(vnl_vector<double> const& x, double *f, vnl_vector<double>* g) = 0;
};

/**
 * Pure affine cost function - parameters are elements of N x N matrix M.
 * Transformation takes place in voxel coordinates - not physical coordinates (for speed)
 */
template <unsigned int VDim, typename TReal = double>
class PureAffineCostFunction : public AbstractAffineCostFunction<VDim, TReal>
{
public:
  typedef AbstractAffineCostFunction<VDim, TReal> Superclass;
  typedef typename Superclass::ParentType ParentType;
  typedef typename Superclass::OFHelperType OFHelperType;
  typedef typename Superclass::LinearTransformType LinearTransformType;

  // Construct the function
  PureAffineCostFunction(GreedyParameters *param, ParentType *parent, int level, OFHelperType *helper);

  // Get the parameters for the specified initial transform
  vnl_vector<double> GetCoefficients(LinearTransformType *tran);

  // Get the transform for the specificed coefficients
  void GetTransform(const vnl_vector<double> &coeff, LinearTransformType *tran);

  // Get the preferred scaling for this function given image dimensions
  virtual vnl_vector<double> GetOptimalParameterScaling(const itk::Size<VDim> &image_dim);

  // Cost function computation
  virtual void compute(vnl_vector<double> const& x, double *f, vnl_vector<double>* g);

protected:
  typedef typename ParentType::ImageType ImageType;
  typedef typename ParentType::ImagePointer ImagePointer;
  typedef typename ParentType::VectorImageType VectorImageType;
  typedef typename ParentType::VectorImagePointer VectorImagePointer;

  // Data needed to compute the cost function
  GreedyParameters *m_Param;
  OFHelperType *m_OFHelper;
  GreedyApproach<VDim, TReal> *m_Parent;
  bool m_Allocated;
  int m_Level;

  // Storage for the gradient of the similarity map
  VectorImagePointer m_Phi, m_GradMetric, m_GradMask;
  ImagePointer m_Metric, m_Mask;

  // Last set of coefficients evaluated
  vnl_vector<double> last_coeff;
};

/**
 * Physical space affine cost function - parameters are elements of affine transform in
 * physical RAS space.
 */
template <unsigned int VDim, typename TReal = double>
class PhysicalSpaceAffineCostFunction : public AbstractAffineCostFunction<VDim, TReal>
{
public:
  typedef AbstractAffineCostFunction<VDim, TReal> Superclass;
  typedef typename Superclass::ParentType ParentType;
  typedef typename Superclass::OFHelperType OFHelperType;
  typedef typename Superclass::LinearTransformType LinearTransformType;

  PhysicalSpaceAffineCostFunction(GreedyParameters *param, ParentType *parent, int level, OFHelperType *helper);
  virtual vnl_vector<double> GetCoefficients(LinearTransformType *tran);
  virtual void GetTransform(const vnl_vector<double> &coeff, LinearTransformType *tran);
  virtual void compute(vnl_vector<double> const& x, double *f, vnl_vector<double>* g);
  virtual vnl_vector<double> GetOptimalParameterScaling(const itk::Size<VDim> &image_dim);

  void map_phys_to_vox(const vnl_vector<double> &x_phys, vnl_vector<double> &x_vox);

protected:
  PureAffineCostFunction<VDim, TReal> m_PureFunction;

  // Voxel to physical transforms for fixed, moving image
  typedef vnl_matrix_fixed<double, VDim, VDim> Mat;
  typedef vnl_vector_fixed<double, VDim> Vec;

  Mat Q_fix, Q_mov, Q_fix_inv, Q_mov_inv;
  Vec b_fix, b_mov, b_fix_inv, b_mov_inv;

  vnl_matrix<double> J_phys_vox;
};

/** Abstract scaling cost function - wraps around another cost function and provides scaling */
template <unsigned int VDim, typename TReal = double>
class ScalingCostFunction : public AbstractAffineCostFunction<VDim, TReal>
{
public:
  typedef AbstractAffineCostFunction<VDim, TReal> Superclass;
  typedef typename Superclass::ParentType ParentType;
  typedef typename Superclass::OFHelperType OFHelperType;
  typedef typename Superclass::LinearTransformType LinearTransformType;

  // Construct the function
  ScalingCostFunction(Superclass *pure_function, const vnl_vector<double> &scaling)
    : Superclass(pure_function->get_number_of_unknowns()),
      m_PureFunction(pure_function), m_Scaling(scaling) {}

  // Get the parameters for the specified initial transform
  vnl_vector<double> GetCoefficients(LinearTransformType *tran);

  // Get the transform for the specificed coefficients
  void GetTransform(const vnl_vector<double> &coeff, LinearTransformType *tran);

  // Cost function computation
  virtual void compute(vnl_vector<double> const& x, double *f, vnl_vector<double>* g);

  const vnl_vector<double> &GetScaling() { return m_Scaling; }

protected:

  // Data needed to compute the cost function
  Superclass *m_PureFunction;
  vnl_vector<double> m_Scaling;
};

/** Cost function for rigid registration */
template <unsigned int VDim, typename TReal = double>
class RigidCostFunction : public AbstractAffineCostFunction<VDim, TReal>
{
public:
  typedef AbstractAffineCostFunction<VDim, TReal> Superclass;
  typedef typename Superclass::ParentType ParentType;
  typedef typename Superclass::OFHelperType OFHelperType;
  typedef typename Superclass::LinearTransformType LinearTransformType;

  typedef vnl_vector_fixed<double, VDim> Vec;
  typedef vnl_matrix_fixed<double, VDim, VDim> Mat;

  RigidCostFunction(GreedyParameters *param, ParentType *parent, int level, OFHelperType *helper);
  vnl_vector<double> GetCoefficients(LinearTransformType *tran);
  void GetTransform(const vnl_vector<double> &coeff, LinearTransformType *tran);
  virtual void compute(vnl_vector<double> const& x, double *f, vnl_vector<double>* g);

  // Get the preferred scaling for this function given image dimensions
  virtual vnl_vector<double> GetOptimalParameterScaling(const itk::Size<VDim> &image_dim);

  // Generate a random rotation matrix with rotation angle alpha (radians)
  static Mat GetRandomRotation(vnl_random &randy, double alpha);

protected:

  static Mat GetRotationMatrix(const Vec &q);
  static Vec GetAxisAngle(const Mat &R);

  // We wrap around a physical space affine function, since rigid in physical space is not
  // the same as rigid in voxel space
  PhysicalSpaceAffineCostFunction<VDim, TReal> m_AffineFn;

  // Flip matrix -- allows rigid registration with flips. If the input matrix has a flip,
  // that flip is maintained during registration
  Mat flip;

};


/** Specialized function for 2D rigid registration */
template <typename TReal>
class RigidCostFunction<2, TReal> : public AbstractAffineCostFunction<2, TReal>
{
public:
  typedef AbstractAffineCostFunction<2, TReal> Superclass;
  typedef typename Superclass::ParentType ParentType;
  typedef typename Superclass::OFHelperType OFHelperType;
  typedef typename Superclass::LinearTransformType LinearTransformType;

  const static int VDim = 2;
  typedef vnl_vector_fixed<double, VDim> Vec;
  typedef vnl_matrix_fixed<double, VDim, VDim> Mat;


  RigidCostFunction(GreedyParameters *param, ParentType *parent, int level, OFHelperType *helper);
  vnl_vector<double> GetCoefficients(LinearTransformType *tran);
  void GetTransform(const vnl_vector<double> &coeff, LinearTransformType *tran);
  virtual void compute(vnl_vector<double> const& x, double *f, vnl_vector<double>* g);

  // Get the preferred scaling for this function given image dimensions
  virtual vnl_vector<double> GetOptimalParameterScaling(const itk::Size<VDim> &image_dim);

  // Generate a random rotation matrix with rotation angle alpha (radians)
  static Mat GetRandomRotation(vnl_random &randy, double alpha);

protected:

  static Mat GetRotationMatrix(double theta);
  static double GetRotationAngle(const Mat &R);

  // We wrap around a physical space affine function, since rigid in physical space is not
  // the same as rigid in voxel space
  PhysicalSpaceAffineCostFunction<VDim, TReal> m_AffineFn;

  // Flip matrix -- allows rigid registration with flips. If the input matrix has a flip,
  // that flip is maintained during registration
  Mat flip;

};



#endif // AFFINECOSTFUNCTIONS_H
