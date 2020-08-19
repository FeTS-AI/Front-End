/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: SimpleWarpImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2009-10-29 11:19:10 $
  Version:   $Revision: 1.34 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __SimpleWarpImageFilter_txx
#define __SimpleWarpImageFilter_txx
#include "SimpleWarpImageFilter.h"

#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkNumericTraits.h"
#include "itkProgressReporter.h"
#include "itkContinuousIndex.h"
#include "vnl/vnl_math.h"
namespace itk
{

/**
 * Default constructor.
 */
template <class TInputImage,class TOutputImage,class TDeformationField, class TFloat>
SimpleWarpImageFilter<TInputImage,TOutputImage,TDeformationField,TFloat>
::SimpleWarpImageFilter()
{
  // Setup the number of required inputs
  this->SetNumberOfRequiredInputs( 2 );  
  
  // Setup default values
  m_Interpolator = NULL;
  m_EdgePaddingValue = NumericTraits<PixelType>::Zero;
}

/**
 * Standard PrintSelf method.
 */
template <class TInputImage,class TOutputImage,class TDeformationField, class TFloat>
void
SimpleWarpImageFilter<TInputImage,TOutputImage,TDeformationField,TFloat>
::PrintSelf(std::ostream& os, Indent indent) const
{

  Superclass::PrintSelf(os, indent);
  os << indent << "EdgePaddingValue: "
     << static_cast<typename NumericTraits<PixelType>::PrintType>(m_EdgePaddingValue)
     << std::endl;
}

/**
 * Set deformation field as Inputs[1] for this ProcessObject.
 *
 */
template <class TInputImage,class TOutputImage,class TDeformationField, class TFloat>
void
SimpleWarpImageFilter<TInputImage,TOutputImage,TDeformationField,TFloat>
::SetDeformationField(
  const DeformationFieldType * field )
{
  // const cast is needed because the pipeline is not const-correct.
  DeformationFieldType * input =  
       const_cast< DeformationFieldType * >( field );
  this->ProcessObject::SetNthInput( 1, input );
}


/**
 * Return a pointer to the deformation field.
 */
template <class TInputImage,class TOutputImage,class TDeformationField, class TFloat>
typename SimpleWarpImageFilter<TInputImage,TOutputImage,TDeformationField,TFloat>
::DeformationFieldType *
SimpleWarpImageFilter<TInputImage,TOutputImage,TDeformationField,TFloat>
::GetDeformationField(void)
{
  return static_cast<DeformationFieldType *>
    ( this->ProcessObject::GetInput( 1 ));
}


/**
 * Setup state of filter before multi-threading.
 * InterpolatorType::SetInputImage is not thread-safe and hence
 * has to be setup before ThreadedGenerateData
 */
template <class TInputImage,class TOutputImage,class TDeformationField, class TFloat>
void
SimpleWarpImageFilter<TInputImage,TOutputImage,TDeformationField,TFloat>
::BeforeThreadedGenerateData()
{
  if( !m_Interpolator )
    {
    itkExceptionMacro(<< "Interpolator not set");
    }
  DeformationFieldPointer fieldPtr = this->GetDeformationField();

  // Connect input image to interpolator
  m_Interpolator->SetInputImage( this->GetInput() );
  typename DeformationFieldType::RegionType defRegion = 
    fieldPtr->GetLargestPossibleRegion();
  typename OutputImageType::RegionType outRegion =
    this->GetOutput()->GetLargestPossibleRegion();
}

/**
 * Setup state of filter after multi-threading.
 */
template <class TInputImage,class TOutputImage,class TDeformationField, class TFloat>
void
SimpleWarpImageFilter<TInputImage,TOutputImage,TDeformationField,TFloat>
::AfterThreadedGenerateData()
{
  // Disconnect input image from interpolator
  m_Interpolator->SetInputImage( NULL );
}

/**
 * Compute the output for the region specified by outputRegionForThread.
 */
template <class TInputImage,class TOutputImage,class TDeformationField, class TFloat>
void
SimpleWarpImageFilter<TInputImage,TOutputImage,TDeformationField,TFloat>
::ThreadedGenerateData(
  const OutputImageRegionType& outputRegionForThread,
  ThreadIdType threadId )
{
  InputImageConstPointer inputPtr = this->GetInput();
  OutputImagePointer outputPtr = this->GetOutput();
  DeformationFieldPointer fieldPtr = this->GetDeformationField();

  // iterator for the output image
  ImageRegionIteratorWithIndex<OutputImageType> outputIt(
    outputPtr, outputRegionForThread );
  IndexType index;
  itk::ContinuousIndex<TFloat,ImageDimension> cix;
  DisplacementType displacement;
  
  // iterator for the deformation field
  ImageRegionIterator<DeformationFieldType> 
    fieldIt(fieldPtr, outputRegionForThread );

  while( !outputIt.IsAtEnd() )
    {
    // get the output image index
    index = outputIt.GetIndex();

    // get the required displacement
    displacement = fieldIt.Get();

    // compute the required input image point
    for(unsigned int j = 0; j < ImageDimension; j++ )
      {
      cix[j] = index[j] + m_DeformationScaling * displacement[j];
      }

    // get the interpolated value
    if( m_Interpolator->IsInsideBuffer( cix ) )
      {
      PixelType value = 
        static_cast<PixelType>(m_Interpolator->EvaluateAtContinuousIndex( cix ) );
      outputIt.Set( value );
      }
    else
      {
      outputIt.Set( m_EdgePaddingValue );
      }   
    ++outputIt;
    ++fieldIt; 
    }
}


template <class TInputImage,class TOutputImage,class TDeformationField, class TFloat>
void
SimpleWarpImageFilter<TInputImage,TOutputImage,TDeformationField,TFloat>
::GenerateInputRequestedRegion()
{

  // call the superclass's implementation
  Superclass::GenerateInputRequestedRegion();

  // request the largest possible region for the input image
  InputImagePointer inputPtr = 
    const_cast< InputImageType * >( this->GetInput() );

  if( inputPtr )
    {
    inputPtr->SetRequestedRegionToLargestPossibleRegion();
    }

  // just propagate up the output requested region for the 
  // deformation field.
  DeformationFieldPointer fieldPtr = this->GetDeformationField();
  OutputImagePointer outputPtr = this->GetOutput();
  if(fieldPtr.IsNotNull() )
    {
    fieldPtr->SetRequestedRegion( outputPtr->GetRequestedRegion() );
    if(!fieldPtr->VerifyRequestedRegion())
      {
      fieldPtr->SetRequestedRegion(fieldPtr->GetLargestPossibleRegion());
      }
    }
}


template <class TInputImage,class TOutputImage,class TDeformationField, class TFloat>
void
SimpleWarpImageFilter<TInputImage,TOutputImage,TDeformationField,TFloat>
::GenerateOutputInformation()
{
  // call the superclass's implementation of this method
  Superclass::GenerateOutputInformation();

  OutputImagePointer outputPtr = this->GetOutput();
  DeformationFieldPointer fieldPtr = this->GetDeformationField();
  outputPtr->SetSpacing( fieldPtr->GetSpacing() );
  outputPtr->SetOrigin( fieldPtr->GetOrigin() );
  outputPtr->SetDirection( fieldPtr->GetDirection() );
  outputPtr->SetLargestPossibleRegion( fieldPtr->
                                       GetLargestPossibleRegion() );
}


} // end namespace itk

#endif
