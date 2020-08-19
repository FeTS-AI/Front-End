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
#include "AffineCostFunctions.h"
#include "MultiImageRegistrationHelper.h"
#include "AffineTransformUtilities.h"

template <unsigned int VDim, typename TReal>
PureAffineCostFunction<VDim, TReal>
::PureAffineCostFunction(GreedyParameters *param, ParentType *parent, int level, OFHelperType *helper)
  : Superclass(VDim * (VDim + 1))
{
  // Store the data
  m_Param = param;
  m_OFHelper = helper;
  m_Level = level;
  m_Parent = parent;

  // Allocate the working images, but do not allocate. We will allocate on demand because
  // these affine cost functions may be created without needing to do any computation
  m_Allocated = false;

  m_Phi = VectorImageType::New();
  m_Phi->CopyInformation(helper->GetReferenceSpace(level));
  m_Phi->SetRegions(helper->GetReferenceSpace(level)->GetBufferedRegion());

  m_GradMetric = VectorImageType::New();
  m_GradMetric->CopyInformation(helper->GetReferenceSpace(level));
  m_GradMetric->SetRegions(helper->GetReferenceSpace(level)->GetBufferedRegion());

  m_GradMask = VectorImageType::New();
  m_GradMask->CopyInformation(helper->GetReferenceSpace(level));
  m_GradMask->SetRegions(helper->GetReferenceSpace(level)->GetBufferedRegion());

  m_Metric = ImageType::New();
  m_Metric->CopyInformation(helper->GetReferenceSpace(level));
  m_Metric->SetRegions(helper->GetReferenceSpace(level)->GetBufferedRegion());

  m_Mask = ImageType::New();
  m_Mask->CopyInformation(helper->GetReferenceSpace(level));
  m_Mask->SetRegions(helper->GetReferenceSpace(level)->GetBufferedRegion());
}


template <unsigned int VDim, typename TReal>
void
PureAffineCostFunction<VDim, TReal>
::compute(const vnl_vector<double> &x, double *f, vnl_vector<double> *g)
{
  // Form a matrix/vector from x
  typename LinearTransformType::Pointer tran = LinearTransformType::New();

  // Set the components of the transform
  unflatten_affine_transform(x.data_block(), tran.GetPointer());

  // Allocate a vector to hold the per-component metric values
  vnl_vector<double> comp_metric;

  // Allocate the memory if needed
  if(!m_Allocated)
    {
    m_Phi->Allocate();
    m_GradMetric->Allocate();
    m_GradMask->Allocate();
    m_Metric->Allocate();
    m_Mask->Allocate();
    m_Allocated = true;
    }

  // Compute the gradient
  double val = 0.0;

  // The scaling of the metric. For some metrics, we need to change sign (to minimize) and also
  // it is more readable if it is scaled by some large factor
  double metric_scale =
      (m_Param->metric == GreedyParameters::NCC
       || m_Param->metric == GreedyParameters::MI
       || m_Param->metric == GreedyParameters::NMI)
      ? -10000.0 : 1.0;

  // The output metric report
  MultiComponentMetricReport out_metric;

  // Gradient output
  typename LinearTransformType::Pointer grad;
  if(g)
    grad = LinearTransformType::New();

  // Perform actual metric computation
  if(m_Param->metric == GreedyParameters::SSD)
    {
    m_OFHelper->ComputeAffineMSDMatchAndGradient(
          m_Level, tran, m_Metric, m_Mask, m_GradMetric, m_GradMask, m_Phi, out_metric, grad);

    }
  else if(m_Param->metric == GreedyParameters::NCC)
    {
    m_OFHelper->ComputeAffineNCCMatchAndGradient(
          m_Level, tran, array_caster<VDim>::to_itkSize(m_Param->metric_radius),
          m_Metric, m_Mask, m_GradMetric, m_GradMask, m_Phi, out_metric, grad);
    }
  else if(m_Param->metric == GreedyParameters::MI || m_Param->metric == GreedyParameters::NMI)
    {
    m_OFHelper->ComputeAffineMIMatchAndGradient(
          m_Level, m_Param->metric == GreedyParameters::NMI,
          tran, m_Metric, m_Mask, m_GradMetric, m_GradMask, m_Phi, out_metric, grad);
    }

  // Handle the gradient
  if(g)
    {
    flatten_affine_transform(grad.GetPointer(), g->data_block());
    (*g) *= metric_scale;
    }

  // Scale the output metric
  out_metric.Scale(metric_scale);

  // Report the output values
  if(f)
    *f = out_metric.TotalMetric;

  // Has the metric improved?
  if(m_Parent->GetMetricLog().size())
    {
    const std::vector<MultiComponentMetricReport> &log = m_Parent->GetMetricLog().back();
    if(log.size() == 0 || log.back().TotalMetric > out_metric.TotalMetric)
      {
      // Record the metric value
      m_Parent->RecordMetricValue(out_metric);

      // Write out the current iteration transform
      if(m_Param->output_intermediate.length())
        {
        vnl_matrix<double> Q_physical = ParentType::MapAffineToPhysicalRASSpace(*m_OFHelper, m_Level, tran);
        m_Parent->WriteAffineMatrixViaCache(m_Param->output_intermediate, Q_physical);
        }
      }
    }

}

template <unsigned int VDim, typename TReal>
vnl_vector<double>
PureAffineCostFunction<VDim, TReal>
::GetCoefficients(LinearTransformType *tran)
{
  vnl_vector<double> x_true(this->get_number_of_unknowns());
  flatten_affine_transform(tran, x_true.data_block());
  return x_true;
}

template <unsigned int VDim, typename TReal>
void
PureAffineCostFunction<VDim, TReal>
::GetTransform(const vnl_vector<double> &coeff, LinearTransformType *tran)
{
  unflatten_affine_transform(coeff.data_block(), tran);
}

template <unsigned int VDim, typename TReal>
vnl_vector<double>
PureAffineCostFunction<VDim, TReal>
::GetOptimalParameterScaling(const itk::Size<VDim> &image_dim)
{
  // Initialize the scaling vector
  vnl_vector<double> scaling(this->get_number_of_unknowns());

  // Set the scaling of the parameters based on image dimensions. This makes it
  // possible to set tolerances in units of voxels. The order of change in the
  // parameters is comparable to the displacement of any point inside the image
  typename LinearTransformType::MatrixType matrix;
  typename LinearTransformType::OffsetType offset;

  for(int i = 0; i < VDim; i++)
    {
    offset[i] = 1.0;
    for(int j = 0; j < VDim; j++)
      matrix(i, j) = image_dim[j];
    }

  typename LinearTransformType::Pointer transform = LinearTransformType::New();
  transform->SetMatrix(matrix);
  transform->SetOffset(offset);
  flatten_affine_transform(transform.GetPointer(), scaling.data_block());

  return scaling;
}

/**
 * PHYSICAL SPACE COST FUNCTION - WRAPS AROUND AFFINE
 */
template <unsigned int VDim, typename TReal>
PhysicalSpaceAffineCostFunction<VDim, TReal>
::PhysicalSpaceAffineCostFunction(GreedyParameters *param, ParentType *parent, int level, OFHelperType *helper)
  : Superclass(VDim * (VDim + 1)), m_PureFunction(param, parent, level, helper)
{
  // The rigid transformation must be rigid in physical space, not in voxel space
  // So in the constructor, we must compute the mappings from the two spaces
  GetVoxelSpaceToNiftiSpaceTransform(helper->GetReferenceSpace(level), Q_fix, b_fix);
  GetVoxelSpaceToNiftiSpaceTransform(helper->GetMovingReferenceSpace(level), Q_mov, b_mov);

  // Compute the inverse transformations
  Q_fix_inv = vnl_matrix_inverse<double>(Q_fix);
  b_fix_inv = - Q_fix_inv * b_fix;

  Q_mov_inv = vnl_matrix_inverse<double>(Q_mov);
  b_mov_inv = - Q_mov_inv * b_mov;

  // Take advantage of the fact that the transformation is linear in A and b to compute
  // the Jacobian of the transformation ahead of time, and "lazily", using finite differences
  int n = VDim * (VDim + 1);
  J_phys_vox.set_size(n, n);
  vnl_vector<double> x_phys(n, 0), x_vox_0(n), x_vox(n);

  // Voxel parameter vector corresponding to zero transform
  this->map_phys_to_vox(x_phys, x_vox_0);

  // Compute each column of the jacobian
  for(int i = 0; i < n; i++)
    {
    x_phys.fill(0);
    x_phys[i] = 1;
    this->map_phys_to_vox(x_phys, x_vox);
    J_phys_vox.set_column(i, x_vox - x_vox_0);
    }


}

template <unsigned int VDim, typename TReal>
void
PhysicalSpaceAffineCostFunction<VDim, TReal>
::map_phys_to_vox(const vnl_vector<double> &x_phys, vnl_vector<double> &x_vox)
{
  Mat A_phys;
  Vec b_phys;

  // unflatten the input parameters into A and b
  unflatten_affine_transform(x_phys.data_block(), A_phys, b_phys);

  // convert into voxel-space affine transform
  Mat A_vox = Q_mov_inv * A_phys * Q_fix;
  Vec b_vox = Q_mov_inv * (A_phys * b_fix + b_phys) + b_mov_inv;

  // Flatten back
  x_vox.set_size(m_PureFunction.get_number_of_unknowns());
  flatten_affine_transform(A_vox, b_vox, x_vox.data_block());
}


template <unsigned int VDim, typename TReal>
void
PhysicalSpaceAffineCostFunction<VDim, TReal>
::compute(const vnl_vector<double> &x, double *f, vnl_vector<double> *g)
{
  // Map to voxel space
  vnl_vector<double> x_vox(m_PureFunction.get_number_of_unknowns());
  this->map_phys_to_vox(x, x_vox);

  // Do we need the gradient?
  if(g)
    {
    // Compute the function and gradient wrt voxel parameters
    vnl_vector<double> g_vox(m_PureFunction.get_number_of_unknowns());
    m_PureFunction.compute(x_vox, f, &g_vox);

    // Transform voxel-space gradient into physical-space gradient
    *g = J_phys_vox.transpose() * g_vox;
    }
  else
    {
    // Just compute the function
    m_PureFunction.compute(x_vox, f, NULL);
    }
}

template <unsigned int VDim, typename TReal>
vnl_vector<double>
PhysicalSpaceAffineCostFunction<VDim, TReal>
::GetOptimalParameterScaling(const itk::Size<VDim> &image_dim)
{
  // TODO: work out scaling for this
  return m_PureFunction.GetOptimalParameterScaling(image_dim);
}

template <unsigned int VDim, typename TReal>
vnl_vector<double>
PhysicalSpaceAffineCostFunction<VDim, TReal>
::GetCoefficients(LinearTransformType *tran)
{
  // The input transform is in voxel space, we must return parameters in physical space
  Mat A_vox, A_phys;
  Vec b_vox, b_phys;

  itk_matrix_to_vnl_matrix(tran->GetMatrix(), A_vox);
  itk_vector_to_vnl_vector(tran->GetOffset(), b_vox);

  // convert into physical-space affine transform
  A_phys = Q_mov * A_vox * Q_fix_inv;
  b_phys = Q_mov * (b_vox - b_mov_inv) - A_phys * b_fix;

  // Flatten
  vnl_vector<double> x(m_PureFunction.get_number_of_unknowns());
  flatten_affine_transform(A_phys, b_phys, x.data_block());

  return x;
}

template <unsigned int VDim, typename TReal>
void
PhysicalSpaceAffineCostFunction<VDim, TReal>
::GetTransform(const vnl_vector<double> &x, LinearTransformType *tran)
{
  // Get voxel-space tranform corresponding to the parameters x
  vnl_vector<double> x_vox(m_PureFunction.get_number_of_unknowns());
  this->map_phys_to_vox(x, x_vox);

  // Unflatten into a transform
  unflatten_affine_transform(x_vox.data_block(), tran);
}


/**
 * SCALING COST FUNCTION - WRAPS AROUND AFFINE
 */
template <unsigned int VDim, typename TReal>
void
ScalingCostFunction<VDim, TReal>
::compute(const vnl_vector<double> &x, double *f, vnl_vector<double> *g)
{
  // Scale the parameters so they are in unscaled units
  vnl_vector<double> x_scaled = element_quotient(x, m_Scaling);

  // Call the wrapped method
  if(g)
    {
    vnl_vector<double> g_scaled(x_scaled.size());
    m_PureFunction->compute(x_scaled, f, &g_scaled);
    *g = element_quotient(g_scaled, m_Scaling);
    }
  else
    {
    m_PureFunction->compute(x_scaled, f, g);
    }
}

// Get the parameters for the specified initial transform
template <unsigned int VDim, typename TReal>
vnl_vector<double>
ScalingCostFunction<VDim, TReal>
::GetCoefficients(LinearTransformType *tran)
{
  vnl_vector<double> x_true = m_PureFunction->GetCoefficients(tran);
  return element_product(x_true, m_Scaling);
}

// Get the transform for the specificed coefficients
template <unsigned int VDim, typename TReal>
void
ScalingCostFunction<VDim, TReal>
::GetTransform(const vnl_vector<double> &coeff, LinearTransformType *tran)
{
  vnl_vector<double> x_true = element_quotient(coeff, m_Scaling);
  m_PureFunction->GetTransform(x_true, tran);
}



/**
 * RIGID COST FUNCTION - WRAPS AROUND AFFINE
 */
template <unsigned int VDim, typename TReal>
RigidCostFunction<VDim, TReal>
::RigidCostFunction(GreedyParameters *param, ParentType *parent, int level, OFHelperType *helper)
  : Superclass(VDim * 2), m_AffineFn(param, parent, level, helper)
{
  // Store the flipped status of the matrix
  this->flip.set_identity();
}

template <unsigned int VDim, typename TReal>
void
RigidCostFunction<VDim, TReal>
::compute(const vnl_vector<double> &x, double *f, vnl_vector<double> *g)
{
  // Place parameters into q and b
  Vec q, b;
  q[0] = x[0]; q[1] = x[1]; q[2] = x[2];
  b[0] = x[3]; b[1] = x[4]; b[2] = x[5];

  // Compute theta
  double theta = q.magnitude();

  // Predefine the rotation matrix
  Mat R; R.set_identity();

  // Create the Q matrix
  Mat Qmat; Qmat.fill(0.0);
  Qmat(0,1) = -q[2]; Qmat(1,0) =  q[2];
  Qmat(0,2) =  q[1]; Qmat(2,0) = -q[1];
  Qmat(1,2) = -q[0]; Qmat(2,1) =  q[0];

  // Compute the square of the matrix
  Mat QQ = vnl_matrix_fixed_mat_mat_mult(Qmat, Qmat);

  // A small epsilon for which a better approximation is R = I + Q
  double eps = 1.0e-4;
  double a1, a2;

  // When theta = 0, rotation is identity
  if(theta > eps)
    {
    // Compute the constant terms in the Rodriguez formula
    a1 = sin(theta) / theta;
    a2 = (1 - cos(theta)) / (theta * theta);

    // Compute the rotation matrix
    R += a1 * Qmat + a2 * QQ;
    }
  else
    {
    R += Qmat;
    }

  // Now we have a rotation and a translation, convert to parameters for the affine function
  vnl_vector<double> x_affine(m_AffineFn.get_number_of_unknowns());
  flatten_affine_transform(this->flip * R, b, x_affine.data_block());

  // Split depending on whether there is gradient to compute
  if(g)
    {
    // Create a vector to store the affine gradient
    vnl_vector<double> g_affine(m_AffineFn.get_number_of_unknowns());
    m_AffineFn.compute(x_affine, f, &g_affine);

    // Compute the matrices d_Qmat
    Mat d_Qmat[3], d_R[3];
    d_Qmat[0].fill(0); d_Qmat[0](1,2) = -1; d_Qmat[0](2,1) =  1;
    d_Qmat[1].fill(0); d_Qmat[1](0,2) =  1; d_Qmat[1](2,0) = -1;
    d_Qmat[2].fill(0); d_Qmat[2](0,1) = -1; d_Qmat[2](1,0) =  1;

    // Compute partial derivatives of R wrt q
    if(theta > eps)
      {
      // Compute the scaling factors in the Rodriguez formula
      double d_a1 = (theta * cos(theta) - sin(theta)) / (theta * theta * theta);
      double d_a2 = (theta * sin(theta) + 2 * cos(theta) - 2) /
                    (theta * theta * theta * theta);

      // Loop over the coordinate and compute the derivative of the rotation matrix wrt x
      for(int p = 0; p < 3; p++)
        {
        // Compute the gradient of the rotation with respect to q[p]
        d_R[p] = d_a1 * q[p] * Qmat +
                 a1 * d_Qmat[p] +
                 d_a2 * q[p] * vnl_matrix_fixed_mat_mat_mult(Qmat, Qmat)
                 + a2 * (vnl_matrix_fixed_mat_mat_mult(d_Qmat[p], Qmat) +
                         vnl_matrix_fixed_mat_mat_mult(Qmat, d_Qmat[p]));
        }
      }
    else
      {
      for(int p = 0; p < 3; p++)
        d_R[p] = d_Qmat[p];
      }

    // Create a matrix to hold the jacobian
    vnl_matrix<double> jac(m_AffineFn.get_number_of_unknowns(), 6);
    jac.fill(0.0);

    // Zero vector
    Vec zero_vec; zero_vec.fill(0.0);
    Mat zero_mat; zero_mat.fill(0.0);

    // Fill out the jacobian
    for(int p = 0; p < 3; p++)
      {
      // Fill the corresponding column
      vnl_vector<double> jac_col_q(m_AffineFn.get_number_of_unknowns());
      flatten_affine_transform(this->flip * d_R[p], zero_vec, jac_col_q.data_block());
      jac.set_column(p, jac_col_q);

      // Also set column on the right (wrt translation)
      vnl_vector<double> jac_col_b(m_AffineFn.get_number_of_unknowns());
      Vec ep; ep.fill(0.0); ep[p] = 1;
      flatten_affine_transform(zero_mat, ep, jac_col_b.data_block());
      jac.set_column(p+3, jac_col_b);
      }

    // Multiply the gradient by the jacobian
    *g = jac.transpose() * g_affine;
    }
  else
    {
    m_AffineFn.compute(x_affine, f, NULL);
    }
}

template <unsigned int VDim, typename TReal>
vnl_vector<double>
RigidCostFunction<VDim, TReal>
::GetOptimalParameterScaling(const itk::Size<VDim> &image_dim)
{
  // Initialize the scaling vector
  vnl_vector<double> scaling(this->get_number_of_unknowns());

  // Scaling is harder for rotations. The rotation parameters are in units of
  // radians. We must figure out how many radians are equivalent to a point in
  // the image moving by a single voxel. That actually works out to be 1/dim.

  // So we take the average of the image dimensions and use that as scaling
  double mean_dim = 0;
  for(int i = 0; i < VDim; i++)
    mean_dim += image_dim[i] / VDim;
  scaling[0] = scaling[1] = scaling[2] = mean_dim;
  scaling[3] = scaling[4] = scaling[5] = 1.0;

  return scaling;
}

template <unsigned int VDim, typename TReal>
vnl_vector<double>
RigidCostFunction<VDim, TReal>
::GetCoefficients(LinearTransformType *tran)
{
  // This affine transform is in voxel space. We must first map it into physical
  vnl_vector<double> x_aff_phys = m_AffineFn.GetCoefficients(tran);
  Mat A; Vec b;
  unflatten_affine_transform(x_aff_phys.data_block(), A, b);

  // If the determinant of A is negative, we need to use the flip
  if(vnl_determinant(A) < 0.0)
    {
    this->flip(0,0) = -1.0;
    }
  else
    {
    this->flip(0,0) = 1.0;
    }

  // Compute polar decomposition of the affine matrix
  vnl_svd<double> svd(this->flip * A);
  Mat R = svd.U() * svd.V().transpose();
  Vec q = this->GetAxisAngle(R);

  // Make result
  vnl_vector<double> x(6);
  x[0] = q[0]; x[1] = q[1]; x[2] = q[2];
  x[3] = b[0]; x[4] = b[1]; x[5] = b[2];

  return x;
}

template <unsigned int VDim, typename TReal>
typename RigidCostFunction<VDim, TReal>::Vec
RigidCostFunction<VDim, TReal>
::GetAxisAngle(const Mat &R)
{
  double eps = 1e-4;
  double f_thresh = cos(eps);

  // Compute the matrix logarithm of R
  double f = (vnl_trace(R) - 1) / 2;
  Vec q;
  if(f >= f_thresh)
    {
    q[0] = R(2,1) - R(1,2);
    q[1] = R(0,2) - R(2,0);
    q[2] = R(1,0) - R(0,1);
    q *= 0.5;
    }
  else
    {
    double theta = acos(f);
    double sin_theta = sqrt(1 - f * f);
    q[0] = R(2,1) - R(1,2);
    q[1] = R(0,2) - R(2,0);
    q[2] = R(1,0) - R(0,1);
    q *= theta / (2 * sin_theta);
    }

  return q;
}

template <unsigned int VDim, typename TReal>
typename RigidCostFunction<VDim, TReal>::Mat
RigidCostFunction<VDim, TReal>
::GetRotationMatrix(const Vec &q)
{
  // Compute theta
  double theta = q.magnitude();

  // Predefine the rotation matrix
  Mat R; R.set_identity();

  // Create the Q matrix
  Mat Qmat; Qmat.fill(0.0);
  Qmat(0,1) = -q[2]; Qmat(1,0) =  q[2];
  Qmat(0,2) =  q[1]; Qmat(2,0) = -q[1];
  Qmat(1,2) = -q[0]; Qmat(2,1) =  q[0];

  // Compute the square of the matrix
  Mat QQ = vnl_matrix_fixed_mat_mat_mult(Qmat, Qmat);

  // When theta = 0, rotation is identity
  double eps = 1e-4;

  if(theta > eps)
    {
    // Compute the constant terms in the Rodriguez formula
    double a1 = sin(theta) / theta;
    double a2 = (1 - cos(theta)) / (theta * theta);

    // Compute the rotation matrix
    R += a1 * Qmat + a2 * QQ;
    }
  else
    {
    R += Qmat;
    }

  return R;
}

template <unsigned int VDim, typename TReal>
void
RigidCostFunction<VDim, TReal>
::GetTransform(const vnl_vector<double> &x, LinearTransformType *tran)
{
  // Place parameters into q and b
  Vec q, b;
  q[0] = x[0]; q[1] = x[1]; q[2] = x[2];
  b[0] = x[3]; b[1] = x[4]; b[2] = x[5];

  // Get the rotation matrix
  Mat R = this->GetRotationMatrix(q);

  // This gives us the physical space affine matrices. Flatten and map to voxel space
  vnl_vector<double> x_aff_phys(m_AffineFn.get_number_of_unknowns());
  flatten_affine_transform(this->flip * R, b, x_aff_phys.data_block());
  m_AffineFn.GetTransform(x_aff_phys, tran);
}

template <unsigned int VDim, typename TReal>
typename RigidCostFunction<VDim, TReal>::Mat
RigidCostFunction<VDim, TReal>
::GetRandomRotation(vnl_random &randy, double alpha)
{
  // Generate a random axis of rotation. A triple of Gaussian numbers, normalized to
  // unit length gives a uniform distribution over the sphere
  Vec q_axis;
  for(int d = 0; d < VDim; d++)
    q_axis[d] = randy.normal();
  q_axis.normalize();

  // Generate the axis-angle representation of the rotation
  Vec q = q_axis * alpha;

  // Generate a random rotation using given angles
  return GetRotationMatrix(q);
}


/**
 * 2D RIGID COST FUNCTION - WRAPS AROUND AFFINE
 */
template <typename TReal>
RigidCostFunction<2, TReal>
::RigidCostFunction(GreedyParameters *param, ParentType *parent, int level, OFHelperType *helper)
  : Superclass(3), m_AffineFn(param, parent, level, helper)
{
  // Store the flipped status of the matrix
  this->flip.set_identity();
}

template <typename TReal>
void
RigidCostFunction<2, TReal>
::compute(const vnl_vector<double> &x, double *f, vnl_vector<double> *g)
{
  // Place parameters into theta and b
  double theta = x[0];
  Vec b(x[1], x[2]);

  // Compute the rotation matrix
  Mat R = this->GetRotationMatrix(theta);

  // Now we have a rotation and a translation, convert to parameters for the affine function
  vnl_vector<double> x_affine(m_AffineFn.get_number_of_unknowns());
  flatten_affine_transform(this->flip * R, b, x_affine.data_block());

  // Split depending on whether there is gradient to compute
  if(g)
    {
    // Create a vector to store the affine gradient
    vnl_vector<double> g_affine(m_AffineFn.get_number_of_unknowns());
    m_AffineFn.compute(x_affine, f, &g_affine);

    // Compute the matrices d_Qmat (derivative wrt theta)
    Mat d_R;
    d_R(0,0) = -sin(theta); d_R(0,1) =  cos(theta);
    d_R(1,0) = -cos(theta); d_R(1,1) = -sin(theta);

    // Create a matrix to hold the jacobian
    vnl_matrix<double> jac(m_AffineFn.get_number_of_unknowns(), 3);
    jac.fill(0.0);

    // Zero vector
    Vec zero_vec; zero_vec.fill(0.0);
    Mat zero_mat; zero_mat.fill(0.0);

    // Fill out the rotation column
    vnl_vector<double> jac_col_theta(m_AffineFn.get_number_of_unknowns());
    flatten_affine_transform(this->flip * d_R, zero_vec, jac_col_theta.data_block());
    jac.set_column(0, jac_col_theta);

    // Fill out the translation columns
    for(int p = 0; p < VDim; p++)
      {
      // Fill the corresponding column
      vnl_vector<double> jac_col_b(m_AffineFn.get_number_of_unknowns());
      Vec ep; ep.fill(0.0); ep[p] = 1;
      flatten_affine_transform(zero_mat, ep, jac_col_b.data_block());
      jac.set_column(p+1, jac_col_b);
      }

    // Multiply the gradient by the jacobian
    *g = jac.transpose() * g_affine;
    }
  else
    {
    m_AffineFn.compute(x_affine, f, NULL);
    }
}

template <typename TReal>
vnl_vector<double>
RigidCostFunction<2, TReal>
::GetOptimalParameterScaling(const itk::Size<VDim> &image_dim)
{
  // Initialize the scaling vector
  vnl_vector<double> scaling(this->get_number_of_unknowns());

  // Scaling is harder for rotations. The rotation parameters are in units of
  // radians. We must figure out how many radians are equivalent to a point in
  // the image moving by a single voxel. That actually works out to be 1/dim.

  // So we take the average of the image dimensions and use that as scaling
  double mean_dim = 0;
  for(int i = 0; i < VDim; i++)
    mean_dim += image_dim[i] / VDim;

  scaling[0] = mean_dim;
  scaling[1] = scaling[2] = 1.0;

  return scaling;
}

template <typename TReal>
vnl_vector<double>
RigidCostFunction<2, TReal>
::GetCoefficients(LinearTransformType *tran)
{
  // This affine transform is in voxel space. We must first map it into physical
  vnl_vector<double> x_aff_phys = m_AffineFn.GetCoefficients(tran);
  Mat A; Vec b;
  unflatten_affine_transform(x_aff_phys.data_block(), A, b);

  // If the determinant of A is negative, we need to use the flip
  if(vnl_determinant(A) < 0.0)
    {
    this->flip(0,0) = -1.0;
    }
  else
    {
    this->flip(0,0) = 1.0;
    }

  // Compute polar decomposition of the affine matrix
  vnl_svd<double> svd(this->flip * A);
  Mat R = svd.U() * svd.V().transpose();

  // Use the angle
  double theta = this->GetRotationAngle(R);

  // Make result
  vnl_vector<double> x(3);
  x[0] = theta; x[1] = b[0]; x[2] = b[1];

  return x;
}

template <typename TReal>
typename RigidCostFunction<2, TReal>::Mat
RigidCostFunction<2, TReal>
::GetRandomRotation(vnl_random &, double alpha)
{
  return GetRotationMatrix(alpha);
}

template <typename TReal>
typename RigidCostFunction<2, TReal>::Mat
RigidCostFunction<2, TReal>
::GetRotationMatrix(double theta)
{
  Mat R;
  R(0,0) =  cos(theta); R(0,1) = sin(theta);
  R(1,0) = -sin(theta); R(1,1) = cos(theta);
  return R;
}

template <typename TReal>
double
RigidCostFunction<2, TReal>
::GetRotationAngle(const Mat &R)
{
  return atan2(R(0,1), R(0,0));
}


template <typename TReal>
void
RigidCostFunction<2, TReal>
::GetTransform(const vnl_vector<double> &x, LinearTransformType *tran)
{
  // Place parameters into q and b
  double theta = x[0];
  Vec b(x[1], x[2]);

  // Get the rotation matrix
  Mat R = this->GetRotationMatrix(theta);

  // This gives us the physical space affine matrices. Flatten and map to voxel space
  vnl_vector<double> x_aff_phys(m_AffineFn.get_number_of_unknowns());
  flatten_affine_transform(this->flip * R, b, x_aff_phys.data_block());
  m_AffineFn.GetTransform(x_aff_phys, tran);
}


#define greedy_template_inst(ClassName) \
  template class ClassName<2, float>; \
  template class ClassName<3, float>; \
  template class ClassName<4, float>; \
  template class ClassName<2, double>; \
  template class ClassName<3, double>; \
  template class ClassName<4, double>;

greedy_template_inst(PureAffineCostFunction)
greedy_template_inst(PhysicalSpaceAffineCostFunction)
greedy_template_inst(ScalingCostFunction)
greedy_template_inst(RigidCostFunction)
