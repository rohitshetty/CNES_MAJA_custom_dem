/*
* Copyright (C) 2020 Centre National d'Etudes Spatiales (CNES)
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*    http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*
*/
/************************************************************************************************************
 *                                                                                                          *
 *                                ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo         *
 *                             o                                                                            *
 *                          o                                                                               *
 *                        o                                                                                 *
 *                      o                                                                                   *
 *                     o       ooooooo       ooooooo    o         o      oo                                 *
 *    o              o       o        o     o       o   o         o     o   o                               *
 *      o           o       o          o   o         o  o         o    o      o                             *
 *        o        o       o           o   o         o  o         o    o        o                           *
 *         o      o        o      oooo     o         o  o         o   o           o                         *
 *          o    o          o              o         o  o         o   o           o                         *
 *           o  o            o             o         o  o o      o   o          o                           *
 *            oo              oooooooo    o         o   o  oooooo   o      oooo                             *
 *                                                     o                                                    *
 *                                                     o                                                    *
 *                                                    o                            o                        *
 *                                                    o            o      oooo     o   o      oooo          *
 *                                                   o             o         o    o    o         o          *
 *                                                   o            o       ooo     o   o       ooo           *
 *                                                               o       o       o   o          o           *
 *                                                               ooooo   oooo    o   ooooo  oooo            *
 *                                                                              o                           *
 *                                                                                                          *
 ************************************************************************************************************
 *                                                                                                          *
 * Author: CS Systemes d'Information  (France)                                                              *
 *                                                                                                          *
 ************************************************************************************************************
 * HISTORIQUE                                                                                               *
 *                                                                                                          *
 * VERSION : 4-0-0 : FA : LAIG-FA-MAC-121054-CS : 24 avril 2014 : Creation suite ?? la correction de la      *
 *                                                          rasterisation des angles PDV et solaire pour S2 *
 *                                                                                                          *
 * FIN-HISTORIQUE                                                                                           *
 *                                                                                                          *
 * $Id: vnsRawExpandVectorImageFilter.txx 2656 2011-02-15 10:25:38Z cvallade $
 *                                                                                                          *
 ************************************************************************************************************/
#ifndef __vnsRawExpandImageFilter_h
#define __vnsRawExpandImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"

namespace vns
{

/** \class RawExpandImageFilter
 * \brief Expand the size of an image by an integer factor in each
 * dimension.
 *
 * RawExpandImageFilter increases the size of an image by an integer
 * factor in each dimension using a interpolation method.
 * The output image size in each dimension is given by:
 *
 * OutputSize[j] = InputSize[j] * ExpandFactors[j]
 *
 * The output values are obtained by interpolating the input image.
 * The default interpolation type used is the LinearInterpolateImageFunction.
 * The user can specify a particular interpolation function via
 * SetInterpolator(). Note that the input interpolator must derive
 * from base class InterpolateImageFunction.
 *
 * \warning: The following is valid only when the flag
 * ITK_USE_CENTERED_PIXEL_COORDINATES_CONSISTENTLY is OFF.
 * When the LargestPossibleRegion is requested, the output image will contain
 * padding at the upper edge of each dimension. The width of padding in the
 * i'th dimension is (ExpandFactors[i] - 1). Users can specify the padding
 * value used by setting the EdgePaddingValue.
 *
 * \warning: The following is valid only when the flag
 * ITK_USE_CENTERED_PIXEL_COORDINATES_CONSISTENTLY is ON 
 * The output image will not contain any padding, and therefore the
 * EdgePaddingValue will not be used.
 *
 * This filter will produce an output with different pixel spacing
 * that its input image such that:
 *
 * OutputSpacing[j] = InputSpacing[j] / ExpandFactors[j]
 *
 * The filter is templated over the input image type and the output 
 * image type.
 *
 * This filter is implemented as a multithreaded filter and supports
 * streaming.
 *
 * \warning This filter only works for image with scalar pixel types.
 * For vector images use VectorRawExpandImageFilter.
 *
 * This filter assumes that the input and output image has the same
 * number of dimensions.
 *
 * \sa InterpolateImageFunction
 * \sa LinearInterpolationImageFunction
 * \sa VectorRawExpandImageFilter
 *
 * \ingroup GeometricTransform
 */
template <class TInputImage, class TOutputImage>
class RawExpandImageFilter:
    public itk::ImageToImageFilter<TInputImage,TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef RawExpandImageFilter                             Self;
  typedef itk::ImageToImageFilter<TInputImage,TOutputImage>  Superclass;
  typedef itk::SmartPointer<Self>                            Pointer;
  typedef itk::SmartPointer<const Self>                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);  

  /** Run-time type information (and related methods). */
  itkTypeMacro(RawExpandImageFilter, itk::ImageToImageFilter);

  /** Typedef to describe the output image region type. */
  typedef typename TOutputImage::RegionType OutputImageRegionType;

  /** ImageDimension enumeration. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  /** Inherit some types from superclass. */
  typedef typename Superclass::InputImageType  InputImageType;
  typedef typename Superclass::OutputImageType OutputImageType;
  typedef typename OutputImageType::PixelType  OutputPixelType;
  typedef typename InputImageType::Pointer     InputImagePointer;
  typedef typename OutputImageType::Pointer    OutputImagePointer;

  /** Factor type */
  typedef double                             ExpandFactorsType;

  /** Typedef support for the interpolation function. */
  typedef double                             CoordRepType;
  typedef itk::InterpolateImageFunction<InputImageType,CoordRepType>
                                             InterpolatorType;
  typedef typename InterpolatorType::Pointer InterpolatorPointer;
  typedef itk::LinearInterpolateImageFunction<InputImageType,CoordRepType>
                                             DefaultInterpolatorType;

  /** Set the interpolator function. */
  itkSetObjectMacro( Interpolator, InterpolatorType );

  /** Get a pointer to the interpolator function. */
  itkGetObjectMacro( Interpolator, InterpolatorType );

  /** Set the expand factors. Values are clamped to 
   * a minimum value of 1. Default is 1 for all dimensions. */
  virtual void SetExpandFactors( const ExpandFactorsType factors[] );
  virtual void SetExpandFactors( const ExpandFactorsType factor );

  /** Get the expand factors. */
  virtual const ExpandFactorsType * GetExpandFactors() const
  { return m_ExpandFactors; }

  /** Set the edge padding value. The default is zero. */
  itkSetMacro( EdgePaddingValue, OutputPixelType );

  /** Get the edge padding value */
  itkGetConstMacro( EdgePaddingValue, OutputPixelType );

  /** RawExpandImageFilter produces an image which is a different resolution and
   * with a different pixel spacing than its input image.  As such,
   * RawExpandImageFilter needs to provide an implementation for
   * UpdateOutputInformation() in order to inform the pipeline execution model.
   * The original documentation of this method is below.
   * \sa ProcessObject::GenerateOutputInformaton() */
  virtual void GenerateOutputInformation();

  /** RawExpandImageFilter needs a smaller input requested region than the output
   * requested region.  As such, ShrinkImageFilter needs to provide an 
   * implementation for GenerateInputRequestedRegion() in order to inform 
   * the pipeline execution model.  
   * \sa ProcessObject::GenerateInputRequestedRegion() */
  virtual void GenerateInputRequestedRegion();

protected:
  RawExpandImageFilter();
  virtual ~RawExpandImageFilter() {}
  void PrintSelf(std::ostream& os, itk::Indent indent) const;

  /** RawExpandImageFilter is implemented as a multithreaded filter.  Therefore,
   * this implementation provides a ThreadedGenerateData() routine which
   * is called for each processing thread. The output image data is allocated
   * automatically by the superclass prior to calling ThreadedGenerateData().
   * ThreadedGenerateData can only write to the portion of the output image
   * specified by the parameter "outputRegionForThread"
   *
   * \sa ImageToImageFilter::ThreadedGenerateData(),
   *     ImageToImageFilter::GenerateData() */
  virtual
  void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                            itk::ThreadIdType threadId );
  
  /** This method is used to set the state of the filter before 
   * multi-threading. */
  virtual void BeforeThreadedGenerateData();

private:
  RawExpandImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  ExpandFactorsType      m_ExpandFactors[ImageDimension];
  InterpolatorPointer    m_Interpolator;
  OutputPixelType        m_EdgePaddingValue;
 
};

} // end namespace itk

#ifndef VNS_MANUAL_INSTANTIATION
#include "vnsRawExpandImageFilter.txx"
#endif

#endif
