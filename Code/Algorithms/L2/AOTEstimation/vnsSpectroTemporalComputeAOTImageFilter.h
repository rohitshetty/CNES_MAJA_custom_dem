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
 * VERSION : 3.2.0 : DM : LAIG-DM-MAJA-2737-CNES : 07 juin 2018 : CAMS utilises dans AOTEstimation          *
 * VERSION : 1.0.0 : FA : LAIG-FA-MAC-148399-CS : 26 octobre 2016 : Audit code - Correction qualite         *
 * VERSION : 5-1-0 : FA : LAIG-FA-MAC-145739-CS : 27 juin 2016 : Audit code - Supp de la macro ITK_EXPORT   *
 * VERSION : 5-1-0 : FA : LAIG-FA-MAC-144674-CS : 18 mai 2016 : Correction refactoring,menage, etc...       *
 * VERSION : 5-0-0 : FA : LAIG-FA-MAC-143764-CS : 19 avril 2016 : Mise en forme du code (indentation, etc.) *
 * VERSION : 5-0-0 : DM : LAIG-DM-MAC-982-CNES : 8 janvier 2016 : Modification de l???algorithme de d??tection *
 *                                                                  des ombres de nuages                    *
 * VERSION : 4-0-0 : DM : LAIG-DM-MAC-751-CNES : 17 f??vrier 2014 : S??lection des pixels pour l'inversion    *
 * 																	(Spec v2.2)								*
 * VERSION : 3-2-0 : FA : LAIG-FA-MAC-701-CNES : 21 novembre 2013 : Correction de l'information g??ographique*
 * 									dans les produits. Methode GenerateOutputInformation priv??e.			*
 * VERSION : 2-0-0 : DM : LAIG-DM-MAC-146-CNES : 20 oct. 2011 : Cr??ation de la m??thode spectro temporelle        *
 *
 * VERSION : 1-0-0 : <TypeFT> : <NumFT> : 18 juin 2010 : Creation                                                           
 *                                                                                                          *
 * FIN-HISTORIQUE                                                                                           *
 *                                                                                                          *
 * $Id$
 *                                                                                                          *
 ************************************************************************************************************/
#ifndef vnsSpectroTemporalComputeAOTImageFilter_h
#define vnsSpectroTemporalComputeAOTImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkNeighborhoodIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include "vnsSpectroTemporalAOTCostFunction.h"
#include "vnsRunLevenbergMarquardtOptimization.h"
#include "vnsMultiTemporalAOTMethodFunctions.h"
#include "vnsMultiSpectralAOTMethodFunctions.h"
#include "vnsMacro.h"
#include "vnsUtilities.h"
#include <map>

namespace vns
{
    /** \class  SpectroTemporalComputeAOTImageFilter
     * \brief This class implements .
     *
     *
     *
     * \author CS Systemes d'Information
     *
     * \sa ImageToImageFilter
     *
     * TInputImage  : vector image of double expected
     * TInputMask   : Monochannel image of short expected
     * TInputDate   : Monochannel image of short expected
     * TOutputImage : Monochannel image of double expected
     * TOutputMask  : Monochannel image of short expected
     *
     * \ingroup L2
     *
     */
    template<class TInputImage, class TInputMask, class TInputVectorMask, class TInputDate, class TInputDTM, class TInputLUT,
            class TOutputImage, class TOutputMask>
        class SpectroTemporalComputeAOTImageFilter : public itk::ImageToImageFilter<TInputImage, TOutputImage>
        {
        public:
            /** Standard class typedefs. */
            typedef SpectroTemporalComputeAOTImageFilter Self;
            typedef itk::ImageToImageFilter<TInputImage, TOutputImage> Superclass;
            typedef itk::SmartPointer<Self> Pointer;
            typedef itk::SmartPointer<const Self> ConstPointer;

            /** Type macro */
            itkNewMacro(Self)

            /** Creation through object factory macro */
            itkTypeMacro(SpectroTemporalComputeAOTImageFilter, ImageToImageFilter )

            /** Some convenient typedefs. */
            typedef TInputImage InputImageType;
            typedef typename InputImageType::ConstPointer InputImageConstPointerType;
            typedef typename InputImageType::Pointer InputImagePointerType;
            typedef typename InputImageType::RegionType InputImageRegionType;
            typedef typename InputImageType::PixelType InputImagePixelType;
            typedef typename InputImageType::InternalPixelType InternalPixelType;
            typedef typename InputImageType::IndexType IndexType;
            typedef typename InputImageType::SizeType SizeType;
            typedef typename InputImageType::OffsetType OffsetType;

            /** Input AOT image typedef */
            typedef TOutputImage AOTImageType;
            typedef typename AOTImageType::ConstPointer AOTImageConstPointerType;
            typedef typename AOTImageType::Pointer AOTImagePointerType;
            typedef typename AOTImageType::PixelType AOTPixelType;

            typedef TInputMask InputMaskType;
            typedef typename InputMaskType::Pointer InputMaskPointerType;
            typedef typename InputMaskType::ConstPointer InputMaskConstPointerType;
            typedef typename InputMaskType::PixelType InputMaskPixelType;

            typedef TInputVectorMask InputVectorMaskType;
            typedef typename InputVectorMaskType::Pointer InputVectorMaskPointerType;
            typedef typename InputVectorMaskType::ConstPointer InputVectorMaskConstPointerType;
            typedef typename InputVectorMaskType::PixelType InputVectorMaskPixelType;

            typedef TInputDate InputDateType;
            typedef typename InputDateType::Pointer InputDatePointerType;
            typedef typename InputDateType::ConstPointer InputDateConstPointerType;
            typedef typename InputDateType::PixelType InputDatePixelType;
            typedef typename InputDateType::InternalPixelType InputDateInternalPixelType;
            typedef InputDateInternalPixelType DateType;

            typedef TInputDTM InputDTMType;
            typedef typename InputDTMType::Pointer InputDTMPointerType;
            typedef typename InputDTMType::ConstPointer InputDTMConstPointerType;
            typedef typename InputDTMType::PixelType InputDTMPixelType;
            typedef typename InputDTMType::InternalPixelType InputDTMInternalPixelType;

            typedef TInputLUT LUTType;
            typedef typename LUTType::Pointer LUTPointerType;
            typedef typename LUTType::PixelType LUTPixelType;
            typedef typename LUTType::ParameterValuesType ParameterValuesType;
            typedef std::map<unsigned int, LUTPointerType> LUTMapType;

            typedef TOutputImage OutputImageType;
            typedef typename OutputImageType::Pointer OutputImagePointerType;
            typedef typename OutputImageType::ConstPointer OutputImageConstPointerType;
            typedef typename OutputImageType::RegionType OutputImageRegionType;
            typedef typename OutputImageType::InternalPixelType OutputImagePixelType;
            typedef typename OutputImageType::IndexType OutputImageIndexType;
            typedef typename OutputImageType::SizeType OutputImageSizeType;
            typedef OutputImagePixelType NoDataType;

            typedef TOutputMask OutputMaskType;
            typedef typename OutputMaskType::Pointer OutputMaskPointerType;
            typedef typename OutputMaskType::ConstPointer OutputMaskConstPointer;
            typedef typename OutputMaskType::PixelType OutputMaskPixelType;
            typedef typename OutputMaskType::RegionType OutputMaskRegionType;

            typedef itk::ImageBase<2> ImageBaseType;

            typedef Utilities::ListOfUIntegers ListOfUIntegers;
            typedef Utilities::ListOfDoubles ListOfDoubles;
            typedef std::vector<InputImagePixelType> PixelVectType;
            typedef std::vector<IndexType> IndexVectType;
            typedef std::vector<OffsetType> OffsetVectType;
            typedef std::list<DateType> DateListType;

            typedef typename itk::ImageRegionIteratorWithIndex<OutputImageType> OutputImageRegionIterator;
            typedef typename itk::ImageRegionIteratorWithIndex<OutputMaskType> OutputMaskRegionIterator;

            /** Neighborhood Iterators */
            typedef itk::ConstNeighborhoodIterator<InputImageType> ImageNeighborhoodIteratorType;
            typedef itk::ConstNeighborhoodIterator<InputMaskType> MaskNeighborhoodIteratorType;
            typedef itk::ConstNeighborhoodIterator<InputVectorMaskType> VectorMaskNeighborhoodIteratorType;
            typedef itk::ConstNeighborhoodIterator<InputDateType> DateNeighborhoodIteratorType;
            typedef itk::ConstNeighborhoodIterator<InputDTMType> DTMNeighborhoodIteratorType;

            typedef typename ImageNeighborhoodIteratorType::RadiusType RadiusType;

            /** FaceLists */
            typedef itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType> ImageBoundaryFacesType;
            typedef typename ImageBoundaryFacesType::FaceListType ImageFaceListType;
            typedef typename ImageFaceListType::iterator ImageFaceListIteratorType;
            typedef typename itk::ZeroFluxNeumannBoundaryCondition<InputImageType> ImageZeroFluxNeumannBCType;

            typedef itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputMaskType> MaskBoundaryFacesType;
            typedef typename MaskBoundaryFacesType::FaceListType MaskFaceListType;
            typedef typename MaskFaceListType::iterator MaskFaceListIteratorType;
            typedef typename itk::ZeroFluxNeumannBoundaryCondition<InputMaskType> MaskZeroFluxNeumannBCType;

            typedef itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputVectorMaskType> VectorMaskBoundaryFacesType;
            typedef typename VectorMaskBoundaryFacesType::FaceListType VectorMaskFaceListType;
            typedef typename VectorMaskFaceListType::iterator VectorMaskFaceListIteratorType;
            typedef typename itk::ZeroFluxNeumannBoundaryCondition<InputVectorMaskType> VectorMaskZeroFluxNeumannBCType;

            typedef itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputDateType> DateBoundaryFacesType;
            typedef typename DateBoundaryFacesType::FaceListType DateFaceListType;
            typedef typename DateFaceListType::iterator DateFaceListIteratorType;
            typedef typename itk::ZeroFluxNeumannBoundaryCondition<InputDateType> DateZeroFluxNeumannBCType;

            typedef itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputDTMType> DTMBoundaryFacesType;
            typedef typename DTMBoundaryFacesType::FaceListType DTMFaceListType;
            typedef typename DTMFaceListType::iterator DTMFaceListIteratorType;
            typedef typename itk::ZeroFluxNeumannBoundaryCondition<InputDTMType> DTMZeroFluxNeumannBCType;

            typedef SpectroTemporalAOTCostFunction<InputImageType, InputDTMType, LUTType> AOTCostFunctionType;

            typedef RunLevenbergMarquardtOptimization<InputImageType, InputDTMType, LUTType, OutputImageType, AOTCostFunctionType> RunLMOptimizationType;
            typedef typename RunLMOptimizationType::Pointer RunLMOptimizationPointerType;

            typedef MultiTemporalAOTMethodFunctions<InputImageType, InputMaskType, InputVectorMaskType, InputDateType> MultiTemporalAOTMethodFunctionsType;

            typedef MultiSpectralAOTMethodFunctions<InputImageType, InputMaskType, InputVectorMaskType> MultiSpectralAOTMethodFunctionsType;

            /** Setters & Getters : Overload SetInput method. */
            vnsSetGetInputRawMacro(L2CLDSub, InputMaskType, 1)
            vnsSetGetInputRawMacro(IPSHDSub, InputMaskType, 2)
            vnsSetGetInputRawMacro(IPSAT, InputVectorMaskType, 3)
            vnsSetGetInputRawMacro(L2WASSub, InputMaskType, 4)
            vnsSetGetInputRawMacro(IPEDG, InputMaskType, 5)
            vnsSetGetInputRawMacro(IPSNWSub, InputMaskType, 6)
            vnsSetGetInputRawMacro(L2inNDT, InputMaskType, 7)
            vnsSetGetInputRawMacro(IPTOACSub, InputImageType, 0)
            vnsSetGetInputRawMacro(L2inRTA, InputImageType, 8)
            vnsSetGetInputRawMacro(L2inRTC, InputImageType, 9)
            vnsSetGetInputRawMacro(L2inPXD, InputDateType, 10)
            vnsSetGetInputRawMacro(DTMSub, InputDTMType, 11)
			vnsSetGetInputRawMacro(AOTCamsSub, AOTImageType, 12)

            OutputImageType *
            GetIPAOT();
            OutputImageType *
			GetIPAOTCams();
            OutputMaskType *
            GetOutputMask();

            /** Product julian date accessors */
            itkSetMacro(JDay, double)
            itkGetConstReferenceMacro(JDay, double)

            /** No_data pixel value accessors for output image AOT */
            itkSetMacro(NoData, NoDataType)
            itkGetConstReferenceMacro(NoData, NoDataType)

            /** Minimum radius for the CurrentRegion */
            itkSetMacro(MinSamplingInterval,unsigned int)

            /** Maximum radius for the CurrentRegion */
            itkSetMacro(MaxSamplingInterval,unsigned int)

            /** Increment N distance : step for the sampling interval increase */
            itkSetMacro(IncrementNDistance,unsigned int)

            /** Radius such as the neighborhood size is initialized with :
             * (2*m_Radius+1)*(2*m_Radius+1) */
            itkSetMacro(Radius,unsigned int)

            /** stepAOT * coarse_res = AOT_resolution */
            itkSetMacro(StepAOT,unsigned int)
            itkGetConstMacro(StepAOT,unsigned int)

            /** Maximum NIR variance parameter */
            itkSetMacro(RelativeVarMax,double)

            /** Composite maximum age threshold */
            itkSetMacro(MaxAge,double)

            /** Red channel ID */
            itkSetMacro(RedChannel,unsigned int)

            /** VAR channel ID (11 by default)  for Date D of the product*/
            itkSetMacro(VARBandD,unsigned int)

            /** VAR channel ID (11 by default)  for Date D-1 of the product*/
            itkSetMacro(VARBandDm1,unsigned int)

            /** STDEV channel ID */
            itkSetMacro(STDEVBand,unsigned int)

            /** Std Noise parameter */
            itkSetMacro(StdNoise,double)

            /** Minimum number of pixels parameter */
            itkSetMacro(MinNumberOfPixels,unsigned int)

            /** Standard Deviation threshold parameter */
            itkSetMacro(MinStdDevThreshold,double)

            /** Thresholds used to select pixels int he multi spectral selection */
            itkSetMacro(RedTOAThreshold,double)
            itkSetMacro(NDVIThreshold,double)

            /** Threshold on equation 1 weight coefficient */
            itkSetMacro(MTWeightThreshold,double)

            /** Numerator coefficient used to compute weight of the multi temporal equation */
            itkSetMacro(ACoefficient,double)

            /** Denominator coefficient used to compute weight of the multi temporal equation */
            itkSetMacro(BCoefficient,double)

            /** Weight for inversion equation 2 in the multi temporal method */
            itkSetMacro(Weq2,double)

            /** Vector of all used indexes considering the step StepAOT */
            vnsSetMacro(AllIndexes,IndexVectType)

            /** Low TOA reflectance thresholds */
            vnsSetMacro(LowReflectanceThresholds, ListOfDoubles)

            /** High TOA reflectance thresholds */
            vnsSetMacro(HighReflectanceThresholds, ListOfDoubles)

            /** Levenberg Marquart Optimization parameters */
            itkSetMacro(F_Tolerance,double)
            // Function value tolerance
            itkSetMacro(G_Tolerance,double)
            // Gradient magnitude tolerance
            itkSetMacro(X_Tolerance,double)
            // Search space tolerance
            itkSetMacro(Epsilon_Function,double)
            // Step
            itkSetMacro(Max_Iterations,int)
            // Maximum number of iterations
            itkSetMacro(AOTinitialD,double)
            itkSetMacro(AOTinitialDm1,double)

			//Boolean to activate CAMS
			itkSetMacro(UseCams,bool);
            //error ponderation for CAMS
            itkSetMacro(KPondCAMS,double);
            //Mean Altitude for CAMS
            itkSetMacro(MeanAltitude,double);
            //Height Scale for CAMS computation
            itkSetMacro(HeightScale,double);

            /** Set LUT to the AOTCostFunction object */
            void
            SetDarkestPixelAOT(const double pDarkestPixelAOT)
            {
                m_RunLMOptimization->GetCostFunction()->SetDarkestPixelAOT(pDarkestPixelAOT);
            }

            /** Set LUT to the AOTCostFunction object */
            void
            SetLUT(LUTType * pLUT)
            {
                m_RunLMOptimization->GetCostFunction()->SetLUTD(pLUT);
            }

            /** map of LUT associated to the date LPD */
            void
            SetLUTMap(const LUTMapType& pLUTMap)
            {
                m_LUTMap = pLUTMap;
            }

            void
            SetMTAOTbandsD(const ListOfUIntegers & pAOTbands)
            {
                m_RunLMOptimization->GetCostFunction()->SetAOTbandsD(pAOTbands);
                m_MTAOTbandsD = pAOTbands;
            }
            void
            SetMTAOTbandsReflectanceCompositeDm1(const ListOfUIntegers & pAOTbands)
            {
                m_RunLMOptimization->GetCostFunction()->SetAOTbandsReflectanceCompositeDm1(pAOTbands);
                m_MTAOTbandsReflectanceCompositeDm1 = pAOTbands;
            }
            void
            SetMTAOTbandsLUTCompositeDm1(const ListOfUIntegers & pAOTbands)
            {
                m_RunLMOptimization->GetCostFunction()->SetAOTbandsLUTCompositeDm1(pAOTbands);
                m_MTAOTbandsLUTCompositeDm1 = pAOTbands;
            }

            void
            SetNDVIAOTBand1(const unsigned int pBand)
            {
                m_RunLMOptimization->GetCostFunction()->SetNDVIAOTBand1(pBand);
                m_NDVIAOTBand1 = pBand;
            }

            void
            SetNDVIAOTBand2(const unsigned int pBand)
            {
                m_RunLMOptimization->GetCostFunction()->SetNDVIAOTBand2(pBand);
                m_NDVIAOTBand2 = pBand;
            }

            void
            SetMSAOTbands(const ListOfUIntegers & pAOTbands)
            {
                m_RunLMOptimization->GetCostFunction()->SetMSAOTbands(pAOTbands);
                m_MSAOTbands = pAOTbands;
            }

            void
            SetSlopeMS(const double pSlope)
            {
                m_RunLMOptimization->GetCostFunction()->SetSlopeMS(pSlope);
            }

            void
            SetYInterceptMS(const double pYIntercept)
            {
                m_RunLMOptimization->GetCostFunction()->SetYInterceptMS(pYIntercept);
            }

            void
            SetSlopeNDVI(const double pSlope)
            {
                m_RunLMOptimization->GetCostFunction()->SetSlopeNDVI(pSlope);
            }

            void
            SetYInterceptNDVI(const double pYIntercept)
            {
                m_RunLMOptimization->GetCostFunction()->SetYInterceptNDVI(pYIntercept);
            }

            void
            SetWeqMS(const double pWeqMS)
            {
                m_RunLMOptimization->GetCostFunction()->SetWeqMS(pWeqMS);
            }

            void
            SetWdark(const double pWdark)
            {
                m_RunLMOptimization->GetCostFunction()->SetWdark(pWdark);
            }

            void
            SetWbounds(const double pWbounds)
            {
                m_RunLMOptimization->GetCostFunction()->SetWbounds(pWbounds);
            }

            void
            SetAOTmin(const double pAOTmin)
            {
                m_RunLMOptimization->GetCostFunction()->SetAOTmin(pAOTmin);
            }

            void
            SetPixelList(const PixelVectType & pPixelVect)
            {
                m_RunLMOptimization->GetCostFunction()->SetPixelList(pPixelVect);
            }

        protected:
            /** Constructor */
            SpectroTemporalComputeAOTImageFilter();
            /** Destructor */
            virtual
            ~SpectroTemporalComputeAOTImageFilter();
            /** Use a marging in the requested regions to be able to compute in a
             * changing currentRegion cf. ThreadedGenerateData(). */
            void
            GenerateInputRequestedRegion();
            /** Validate the presence of all ten inputs. If one or more inputs
             * are missing, throw an exception. */
            void
            BeforeThreadedGenerateData();
            /** Multi-thread version GenerateData. */
            void
            ThreadedGenerateData(const InputImageRegionType& outputRegionForThread, itk::ThreadIdType threadId);

            /** Modify output information */
            void
            GenerateOutputInformation();

        private:
            SpectroTemporalComputeAOTImageFilter(const Self&); //purposely not implemented
            void
            operator=(const Self&); //purposely not implemented

            /** Add margings to the requested regions */
            void
            PadByRadiusRequestedRegion(ImageBaseType * pImagePtr);

            /** Call the Levenberg Marquardt least squares inversion */
            OutputImagePixelType
            LeastSquaresInversion(const double pWeq1, const double pWeq2, const unsigned int pLPDconst,
                    const OffsetVectType & pMTExtractedOffsets, const OffsetVectType & pMSExtractedOffsets,
                    ImageNeighborhoodIteratorType & pIPTOACSubIt, ImageNeighborhoodIteratorType & pL2inRTAIt,
                    ImageNeighborhoodIteratorType & pL2inRTCIt, DTMNeighborhoodIteratorType & pDTMSubIt,
					const AOTPixelType &pAotCams,RunLMOptimizationPointerType & pLMOptimizer);

            /** Copy attributes from m_RunLMOptimization to the other optimization objects
             * contained into m_VectorOptim. Used to work with treaded generate data.
             */
            void
            CopyLMOptimizationAttributes();

            void
            CallCopyInputRegionToOutputRegion(OutputImageRegionType & destRegion, const InputImageRegionType & srcRegion);

            void
            CallCopyOutputRegionToInputRegion(InputImageRegionType & destRegion, const OutputImageRegionType & srcRegion);

            /** Julian day of the current product */
            double m_JDay;

            /** No_data value set in the output */
            NoDataType m_NoData;

            /** Minimum sampling interval */
            unsigned int m_MinSamplingInterval;

            /** Maximum sampling interval */
            unsigned int m_MaxSamplingInterval;

            /** Step for the sampling interval increase */
            unsigned int m_IncrementNDistance;

            /** Fixed Radius (spec N_AOT) */
            unsigned int m_Radius;

            /** StepAOT */
            unsigned int m_StepAOT;

            /** Maximum NIR variance parameter */
            double m_RelativeVarMax;

            /** Composite maximum age threshold */
            double m_MaxAge;

            /** Red channel ID (Set in m_MSAOTBands[1]) */
            unsigned int m_RedChannel;

            /** NDVI1_aero band ID (B11 by default) */
            unsigned int m_NDVIAOTBand1;

            /** NDVI2_aero band ID (B11 by default) */
            unsigned int m_NDVIAOTBand2;

            /** VAR channel ID (11 by default) for Date D of the product */
            unsigned int m_VARBandD;

            /** VAR channel ID (11 by default) for Date D-1 of the product  */
            unsigned int m_VARBandDm1;

            /** STDEV channel ID */
            unsigned int m_STDEVBand;

            /** Std Noise parameter */
            double m_StdNoise;

            /** Minimum number of pixels parameter */
            unsigned int m_MinNumberOfPixels;

            /** Standard Deviation threshold parameter */
            double m_MinStdDevThreshold;

            /** Thresholds used to select pixels in the multi spectral selection */
            double m_RedTOAThreshold;
            double m_NDVIThreshold;

            /** Threshold on equation 1 weight coefficient */
            double m_MTWeightThreshold;

            /** Numerator coefficient used to compute weight of the multi temporal equation */
            double m_ACoefficient;

            /** Denominator coefficient used to compute weight of the multi temporal equation */
            double m_BCoefficient;

            /** Weight for inversion equation 2 in the multi temporal method */
            double m_Weq2;

            /** list of the used band for the multi temporal AOT Estimation for the date D */
            ListOfUIntegers m_MTAOTbandsD;

            /** list of the used band for the multi temporal AOT Estimation for the date D-1 for 'reflectance' composite */
            ListOfUIntegers m_MTAOTbandsReflectanceCompositeDm1;

            /** list of the used band for the multi temporal AOT Estimation for the date D-1 for 'lut' composite */
            ListOfUIntegers m_MTAOTbandsLUTCompositeDm1;

            /** list of the used band for the multi spectral AOT Estimation */
            ListOfUIntegers m_MSAOTbands;

            /** Levenberg Marquardt Launcher */
            RunLMOptimizationPointerType m_RunLMOptimization;

            std::vector<RunLMOptimizationPointerType> m_VectorOptim;

            /** Vector of all used indexes considering the step StepAOT*/
            IndexVectType m_AllIndexes;

            /** Low TOA reflectance thresholds */
            ListOfDoubles m_LowReflectanceThresholds;

            /** High TOA reflectance thresholds */
            ListOfDoubles m_HighReflectanceThresholds;

            /** Map of LUT for LPD */
            LUTMapType m_LUTMap;

            double m_F_Tolerance; // Function value tolerance
            double m_G_Tolerance; // Gradient magnitude tolerance
            double m_X_Tolerance; // Search space tolerance
            double m_Epsilon_Function; // Step
            int m_Max_Iterations; // Maximum number of iterations
            double m_AOTinitialD;
            double m_AOTinitialDm1;

            //Boolean to activate CAMS
            bool m_UseCams;
            //error ponderation for CAMS
            double m_KPondCAMS;
            //Mean Altitude for CAMS
            double m_MeanAltitude;
            //Height Scale for CAMS computation
            double m_HeightScale;

        };

} // End namespace vns
#ifndef VNS_MANUAL_INSTANTIATION
#include "vnsSpectroTemporalComputeAOTImageFilter.txx"
#endif

#endif /* vnsSpectroTemporalComputeAOTImageFilter_h */
