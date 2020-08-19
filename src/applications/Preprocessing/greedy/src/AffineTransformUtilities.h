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
#ifndef AFFINETRANSFORMUTILITIES_H
#define AFFINETRANSFORMUTILITIES_H

#include <itkMatrixOffsetTransformBase.h>
#include <itkImageBase.h>

/**
 * Flatten an affine transform to a flat array
 */
template<class TFloat, class TFloatArr, unsigned int VDim>
static void flatten_affine_transform(
    const itk::MatrixOffsetTransformBase<TFloat, VDim, VDim> *transform,
    TFloatArr *flat_array)
{
  int pos = 0;
  for(int i = 0; i < VDim; i++)
    {
    flat_array[pos++] = transform->GetOffset()[i];
    for(int j = 0; j < VDim; j++)
      flat_array[pos++] = transform->GetMatrix()(i,j);
    }
}

template<class TFloat, class TFloatArr, unsigned int VDim>
static void flatten_affine_transform(
    const vnl_matrix_fixed<TFloat, VDim, VDim> &matrix,
    const vnl_vector_fixed<TFloat, VDim> &offset,
    TFloatArr *flat_array)
{
  int pos = 0;
  for(int i = 0; i < VDim; i++)
    {
    flat_array[pos++] = offset[i];
    for(int j = 0; j < VDim; j++)
      flat_array[pos++] = matrix(i,j);
    }
}

/**
 * Unflatten a flat array to an affine transform
 */
template<class TFloat, class TFloatArr, unsigned int VDim>
static void unflatten_affine_transform(
   const TFloatArr *flat_array,
   itk::MatrixOffsetTransformBase<TFloat, VDim, VDim> *transform,
   double scaling = 1.0)
{
  typename itk::MatrixOffsetTransformBase<TFloat, VDim, VDim>::MatrixType matrix;
  typename itk::MatrixOffsetTransformBase<TFloat, VDim, VDim>::OffsetType offset;

  int pos = 0;
  for(int i = 0; i < VDim; i++)
    {
    offset[i] = flat_array[pos++] * scaling;
    for(int j = 0; j < VDim; j++)
      matrix(i, j) = flat_array[pos++] * scaling;
    }

  transform->SetMatrix(matrix);
  transform->SetOffset(offset);
}

template<class TFloat, class TFloatArr, unsigned int VDim>
static void unflatten_affine_transform(
    const TFloatArr *flat_array,
    vnl_matrix_fixed<TFloat, VDim, VDim> &matrix,
    vnl_vector_fixed<TFloat, VDim> &offset,
    double scaling = 1.0)
{
  int pos = 0;
  for(int i = 0; i < VDim; i++)
    {
    offset[i] = flat_array[pos++] * scaling;
    for(int j = 0; j < VDim; j++)
      matrix(i, j) = flat_array[pos++] * scaling;
    }
}

template <class TFloat, unsigned int VDim>
static void set_affine_transform(
    const vnl_matrix_fixed<TFloat, VDim, VDim> &matrix,
    const vnl_vector_fixed<TFloat, VDim> &offset,
    itk::MatrixOffsetTransformBase<TFloat, VDim, VDim> *transform)
{
  typename itk::MatrixOffsetTransformBase<TFloat, VDim, VDim>::MatrixType tmatrix;
  typename itk::MatrixOffsetTransformBase<TFloat, VDim, VDim>::OffsetType toffset;

  for(int i = 0; i < VDim; i++)
    {
    toffset[i] = offset[i];
    for(int j = 0; j < VDim; j++)
      tmatrix(i, j) = matrix(i,j);
    }

  transform->SetMatrix(tmatrix);
  transform->SetOffset(toffset);
}

// Helper function to map from ITK coordiante space to RAS space
template<unsigned int VDim, class TMat, class TVec>
void
GetVoxelSpaceToNiftiSpaceTransform(itk::ImageBase<VDim> *image,
                                   TMat &A,
                                   TVec &b)
{
  // Generate intermediate terms
  typedef typename TMat::element_type TReal;
  vnl_matrix<double> m_dir, m_ras_matrix;
  vnl_diag_matrix<double> m_scale, m_lps_to_ras;
  vnl_vector<double> v_origin, v_ras_offset;

  // Compute the matrix
  m_dir = image->GetDirection().GetVnlMatrix();
  m_scale.set(image->GetSpacing().GetVnlVector());
  m_lps_to_ras.set(vnl_vector<double>(VDim, 1.0));
  m_lps_to_ras[0] = -1;
  m_lps_to_ras[1] = -1;
  A = m_lps_to_ras * m_dir * m_scale;

  // Compute the vector
  v_origin = image->GetOrigin().GetVnlVector();
  b = m_lps_to_ras * v_origin;
}

template <class TITKMatrix, class TVNLMatrix>
void vnl_matrix_to_itk_matrix(
    const TVNLMatrix &vmat,
    TITKMatrix &imat)
{
  for(int r = 0; r < TITKMatrix::RowDimensions; r++)
    for(int c = 0; c < TITKMatrix::ColumnDimensions; c++)
      imat(r,c) = static_cast<typename TITKMatrix::ValueType>(vmat(r,c));
}

template <class TITKVector, class TVNLVector>
void vnl_vector_to_itk_vector(
    const TVNLVector &vvec,
    TITKVector &ivec)
{
  for(int r = 0; r < TITKVector::Dimension; r++)
    ivec[r] = static_cast<typename TITKVector::ValueType>(vvec(r));
}

template <class TITKMatrix, class TVNL>
void itk_matrix_to_vnl_matrix(
    const TITKMatrix &imat,
    vnl_matrix_fixed<TVNL,TITKMatrix::RowDimensions,TITKMatrix::ColumnDimensions>  &vmat)
{
  for(int r = 0; r < TITKMatrix::RowDimensions; r++)
    for(int c = 0; c < TITKMatrix::ColumnDimensions; c++)
      vmat(r,c) = static_cast<TVNL>(imat(r,c));
}

template <class TITKMatrix, class TVNL>
void itk_matrix_to_vnl_matrix(
    const TITKMatrix &imat,
    vnl_matrix<TVNL>  &vmat)
{
  vmat.set_size(TITKMatrix::RowDimensions,TITKMatrix::ColumnDimensions);
  for(int r = 0; r < TITKMatrix::RowDimensions; r++)
    for(int c = 0; c < TITKMatrix::ColumnDimensions; c++)
      vmat(r,c) = static_cast<TVNL>(imat(r,c));
}

template <class TITKVector, class TVNL>
void itk_vector_to_vnl_vector(
    const TITKVector &ivec,
    vnl_vector_fixed<TVNL,TITKVector::Dimension> &vvec)
{
  for(int r = 0; r < TITKVector::Dimension; r++)
    vvec(r) = static_cast<TVNL>(ivec[r]);
}

template <class TITKVector, class TVNL>
void itk_vector_to_vnl_vector(
    const TITKVector &ivec,
    vnl_vector<TVNL> &vvec)
{
  vvec.set_size(TITKVector::Dimension);
  for(int r = 0; r < TITKVector::Dimension; r++)
    vvec(r) = static_cast<TVNL>(ivec[r]);
}

#endif // AFFINETRANSFORMUTILITIES_H
