/*==========================================================================================
 *
 *	Program:	TOMI
 *	Module:		BrainAlignment
 *	Language: C++
 *	Data:			$Date: 2014-03-17 10:47:53$
 *	Version:	$Revision: 0.1$
 *
 *	Copyright Philips Research China. All rights reserved.
 *
 *  This program aligns the brain image into the mni atlas space using affine registration.
 *
 *=========================================================================================*/

#ifndef _BrainAlignmentFilter_cxx_
#define _BrainAlignmentFilter_cxx_

#include "BrainAlignmentFilter.h"

#include <itkCenteredVersorTransformInitializer.h>
#include <itkMattesMutualInformationImageToImageMetric.h>
#include "itkVersorRigid3DTransformOptimizer.h"
#include "itkImageRegistrationMethod.h"
#include "itkMultiResolutionImageRegistrationMethod.h"
#include "itkMultiResolutionPyramidImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkResampleImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkIdentityTransform.h"
#include "itkTimeProbe.h"
#include "itkScaleSkewVersor3DTransform.h"
#include "itkImageMaskSpatialObject.h"

#include "itkDualResampleImageFilter.h"


template<class TImage>
BrainAlignment<TImage>::BrainAlignment()
{
	m_Atlas = 0;
	m_InputImage = 0;
	m_RigidOutputImage = 0;
	m_AffineOutputImage = 0;
	m_AffineCSF = 0;
	m_AffineGM = 0;
	m_AffineWM = 0;
	m_AffineLabel = 0;

	m_NumberOfSamples = 100000;
	m_DownsampleFactor = 2;
}

template<class TImage>
BrainAlignment<TImage>::~BrainAlignment()
{

}

template<class TImage>
typename BrainAlignment<TImage>::ImagePointerType	BrainAlignment<TImage>::DownsampleImage(ImagePointerType image)
{
	typedef itk::ResampleImageFilter< ImageType, ImageType >	DownsampleFilterType;
  typedef itk::LinearInterpolateImageFunction< ImageType, double >    InterpolatorType;
	typedef itk::IdentityTransform<double>	TransformType;
	TransformType::Pointer idTransform = TransformType::New();

	typename DownsampleFilterType::Pointer downsampler = DownsampleFilterType::New();
  typename InterpolatorType::Pointer downsampleInterpolator = InterpolatorType::New();

	typename ImageType::SizeType inputSize = image->GetBufferedRegion().GetSize();
	typename ImageType::SizeType downsampleSize;
	downsampleSize[0] = inputSize[0] / m_DownsampleFactor;
	downsampleSize[1] = inputSize[1] / m_DownsampleFactor;
	downsampleSize[2] = inputSize[2] / m_DownsampleFactor;

	typename ImageType::SpacingType inputSpacing = image->GetSpacing();
	typename ImageType::SpacingType downsampleSpacing;
	downsampleSpacing[0] = inputSpacing[0] / (float)downsampleSize[0] * (float)inputSize[0];
	downsampleSpacing[1] = inputSpacing[1] / (float)downsampleSize[1] * (float)inputSize[1];
	downsampleSpacing[2] = inputSpacing[2] / (float)downsampleSize[2] * (float)inputSize[2];

	downsampler->SetTransform(idTransform);
	downsampler->SetInput( image );
	downsampler->SetInterpolator( downsampleInterpolator );

	downsampler->SetSize(    downsampleSize );
	downsampler->SetOutputOrigin(  image->GetOrigin() );
	downsampler->SetOutputSpacing( downsampleSpacing );
	downsampler->SetOutputDirection( image->GetDirection() );
	downsampler->SetDefaultPixelValue( 0 );
	downsampler->Update();
	
	typename ImageType::Pointer downsampleImage = downsampler->GetOutput();

	return downsampleImage;
	
}

template<class TImage>
void BrainAlignment<TImage>::Update()
{
	// Typical type definitions of registration components
	const unsigned int Dimension = 3;

  typedef itk::VersorRigid3DTransformOptimizer           RigidOptimizerType;

	typedef itk::RegularStepGradientDescentOptimizer       AffineOptimizerType;

  typedef itk::MattesMutualInformationImageToImageMetric< 
                                    InternalImageType, 
                                    InternalImageType >    MetricType;
  

  typedef itk::LinearInterpolateImageFunction< 
                                    InternalImageType,
                                    double          >    InterpolatorType;

  typedef itk::MultiResolutionImageRegistrationMethod< 
                                    InternalImageType, 
                                    InternalImageType >   RegistrationType;

  typedef itk::MultiResolutionPyramidImageFilter<
                                    InternalImageType,
                                    InternalImageType >   FixedImagePyramidType;
  typedef itk::MultiResolutionPyramidImageFilter<
                                    InternalImageType,
                                    InternalImageType >   MovingImagePyramidType;

  typedef itk::CastImageFilter< 
                        ImageType, InternalImageType > FixedCastFilterType;
  typedef itk::CastImageFilter< 
                        ImageType, InternalImageType > MovingCastFilterType;

  typedef itk::CenteredTransformInitializer< RigidTransformType,
                                             ImageType, 
                                             ImageType 
                                                 >  TransformInitializerType;
	typedef itk::ResampleImageFilter< 
                            ImageType, 
                            ImageType >    ResampleFilterType;

  typedef itk::LinearInterpolateImageFunction< 
                                 ImageType, double >  ResampleInterpolatorType;

	typedef itk::ImageFileWriter<ImageType>	ReaderType;
	
	typedef itk::ImageFileWriter<ImageType>	WriterType;

	typedef itk::ImageFileWriter<InternalImageType>	FloatReaderType;

	typedef itk::ResampleImageFilter< 
                            InternalImageType, 
                            InternalImageType >    FloatResampleFilterType;

  typedef itk::LinearInterpolateImageFunction< 
                                 InternalImageType, double >  FloatResampleInterpolatorType;

	// Load the atlas image
	typename ReaderType::Pointer atlasReader = ReaderType::New();
	std::string atlasName = m_AtlasPath + "/mni_icbm152_t2_tal_nlin_sym_09c.nii";
	atlasReader->SetFileName(atlasName);
	try
	{
		atlasReader->Update();
	}
	catch(itk::ExceptionObject & err)
	{
		std::cerr<< err <<std::endl;
	}


	m_Atlas = atlasReader->GetOutput();

	typename ImageType::Pointer downsampledAtlas = DownsampleImage(m_Atlas);

	//ImageType::Pointer downsampledInputImage = DownsampleImage(m_InputImage);


	/* First perform rigid registration to align subject to atlas*/
  MetricType::Pointer         rigidMetric        = MetricType::New();
  RigidOptimizerType::Pointer rigidOptimizer     = RigidOptimizerType::New();
  InterpolatorType::Pointer   rigidInterpolator  = InterpolatorType::New();
  RegistrationType::Pointer   rigidRegistration  = RegistrationType::New();
  
  FixedImagePyramidType::Pointer fixedRigidImagePyramid = 
      FixedImagePyramidType::New();
  MovingImagePyramidType::Pointer movingRigidImagePyramid =
      MovingImagePyramidType::New();

	// Key components of the rigid registration
  rigidRegistration->SetMetric(        rigidMetric        );
  rigidRegistration->SetOptimizer(     rigidOptimizer     );
  rigidRegistration->SetInterpolator(  rigidInterpolator  );
  rigidRegistration->SetFixedImagePyramid( fixedRigidImagePyramid );
  rigidRegistration->SetMovingImagePyramid( movingRigidImagePyramid ); 

  RigidTransformType::Pointer  rigidTransform = RigidTransformType::New();
  rigidRegistration->SetTransform( rigidTransform );

	// Set parameter for the MI metric
	const int numberOfPixels = m_Atlas->GetBufferedRegion().GetNumberOfPixels();
  rigidMetric->SetNumberOfHistogramBins(50);
  rigidMetric->SetNumberOfSpatialSamples(m_NumberOfSamples);


	// Cast the short image type to float image type
  typename FixedCastFilterType::Pointer fixedRigidCaster   = FixedCastFilterType::New();
  typename MovingCastFilterType::Pointer movingRigidCaster = MovingCastFilterType::New();

	fixedRigidCaster->SetInput(  downsampledAtlas );
	movingRigidCaster->SetInput( m_InputImage );

  rigidRegistration->SetFixedImage(    fixedRigidCaster->GetOutput()    );
  rigidRegistration->SetMovingImage(   movingRigidCaster->GetOutput()   );

  fixedRigidCaster->Update();

  rigidRegistration->SetFixedImageRegion( 
       fixedRigidCaster->GetOutput()->GetBufferedRegion() );	

	// Initialize the transform
  //TransformInitializerType::Pointer initializer = 
  //                                        TransformInitializerType::New();

  //initializer->SetTransform(   rigidTransform );
  //initializer->SetFixedImage(  m_Atlas );
  //initializer->SetMovingImage( m_InputImage );
  //initializer->GeometryOn();

  //typedef RigidTransformType::VersorType  VersorType;
  //typedef VersorType::VectorType     VectorType;

  //VersorType     rotation;
  //VectorType     axis;
  //
  //axis[0] = 1;
  //axis[1] = 0;
  //axis[2] = 0;
  //const double angle = 0;
  //rotation.Set(  axis, angle  );
  //rigidTransform->SetRotation( rotation );
  //rigidTransform->SetTranslation( translation );
  //initializer->InitializeTransform();

  rigidRegistration->SetInitialTransformParameters( rigidTransform->GetParameters() );

	// Set optimizer scales and parameters
  typedef RigidOptimizerType::ScalesType       RigidOptimizerScalesType;
  RigidOptimizerScalesType rigidOptimizerScales( rigidTransform->GetNumberOfParameters() );

  rigidOptimizerScales[0] = 1000;
  rigidOptimizerScales[1] = 1000;
  rigidOptimizerScales[2] = 1000;
  rigidOptimizerScales[3] = 1;
  rigidOptimizerScales[4] = 1;
  rigidOptimizerScales[5] = 1;

  rigidOptimizer->SetScales( rigidOptimizerScales );

  rigidOptimizer->SetMaximumStepLength( 2  ); 
  rigidOptimizer->SetMinimumStepLength( 0.005 );

  rigidOptimizer->SetNumberOfIterations( 1500 );
  rigidOptimizer->MinimizeOn();

  rigidRegistration->SetNumberOfLevels( 3 );

	// Start registration optimization
	std::cout<<"Start rigid registration..."<<std::endl;
	itk::TimeProbe clock;
	clock.Start();
  try 
  { 
    rigidRegistration->Update(); 
		std::cout << "Optimizer stop condition: "
          << rigidRegistration->GetOptimizer()->GetStopConditionDescription()
          << std::endl;

  } 
  catch( itk::ExceptionObject & err ) 
  { 
    std::cerr << "ExceptionObject caught !" << std::endl; 
    std::cerr << err << std::endl; 
    exit;
  } 
	clock.Stop();
	std::cout << "Time: " << clock.GetTotal() << std::endl;
  
  RigidOptimizerType::ParametersType rigidParameters = 
                    rigidRegistration->GetLastTransformParameters();

	m_RigidTransform = RigidTransformType::New();
	m_RigidTransform->SetParameters(rigidParameters);
	m_RigidTransform->SetCenter(rigidTransform->GetCenter());

	// Resample the subject image using the resulting rigid transform using downsampled resolution
	typename ResampleFilterType::Pointer rigidResampler = ResampleFilterType::New();
  typename ResampleInterpolatorType::Pointer rigidResampleInterpolator = ResampleInterpolatorType::New();

	rigidResampler->SetTransform( m_RigidTransform );
	rigidResampler->SetInput( m_InputImage );
	rigidResampler->SetInterpolator( rigidResampleInterpolator );

	rigidResampler->SetSize(    downsampledAtlas->GetLargestPossibleRegion().GetSize() );
	rigidResampler->SetOutputOrigin(  downsampledAtlas->GetOrigin() );
	rigidResampler->SetOutputSpacing( downsampledAtlas->GetSpacing() );
	rigidResampler->SetOutputDirection( downsampledAtlas->GetDirection() );
	rigidResampler->SetDefaultPixelValue( 0 );
	rigidResampler->Update();

	m_RigidOutputImage = rigidResampler->GetOutput();

	//// write output for debug
	//typedef itk::ImageFileWriter<ImageType>	WriterType;
	//WriterType::Pointer rigidWriter = WriterType::New();
	//rigidWriter->SetInput(m_RigidOutputImage);
	//rigidWriter->SetFileName("rigidOutput.nii");
	//try
	//{
	//	rigidWriter->Update();
	//}
	//catch(itk::ExceptionObject & err)
	//{
	//	std::cerr<< err << std::endl;
	//}


	/* Second perform affine registration to align atlas to template */


	// Key components of the rigid registration
	MetricType::Pointer						affineMetric        = MetricType::New();
  AffineOptimizerType::Pointer  affineOptimizer     = AffineOptimizerType::New();
  InterpolatorType::Pointer			affineInterpolator  = InterpolatorType::New();
  RegistrationType::Pointer			affineRegistration  = RegistrationType::New();
	AffineTransformType::Pointer  affineTransform			= AffineTransformType::New();

  affineRegistration->SetMetric(        affineMetric        );
  affineRegistration->SetOptimizer(     affineOptimizer     );
  affineRegistration->SetInterpolator(  affineInterpolator  );
  affineRegistration->SetTransform(			affineTransform			);

  affineRegistration->SetInitialTransformParameters( affineTransform->GetParameters() );

	// Set metric parameters
  affineMetric->SetNumberOfHistogramBins(50);
  affineMetric->SetNumberOfSpatialSamples(m_NumberOfSamples);

	// Set optimizer scales and parameters
  typedef AffineOptimizerType::ScalesType						AffineOptimizerScalesType;
  AffineOptimizerScalesType affineOptimizerScales(	affineTransform->GetNumberOfParameters() );
  affineOptimizerScales[0] =  1.0;
  affineOptimizerScales[1] =  1.0;
  affineOptimizerScales[2] =  1.0;
  affineOptimizerScales[3] =  1.0;
  affineOptimizerScales[4] =  1.0;
  affineOptimizerScales[5] =  1.0;
  affineOptimizerScales[6] =  1.0;
  affineOptimizerScales[7] =  1.0;
  affineOptimizerScales[8] =  1.0;
  affineOptimizerScales[9]  =  0.001;
  affineOptimizerScales[10] =  0.001;
  affineOptimizerScales[11] =  0.001;
  affineOptimizer->SetScales( affineOptimizerScales );

  affineOptimizer->SetMaximumStepLength( 1 );
  affineOptimizer->SetMinimumStepLength( 0.005 );
  affineOptimizer->SetNumberOfIterations( 200 );

	affineOptimizer->MinimizeOn();

	// Cast the short image type to float image type
  typename FixedCastFilterType::Pointer	fixedAffineCaster  = FixedCastFilterType::New();
  typename MovingCastFilterType::Pointer movingAffineCaster = MovingCastFilterType::New();

	fixedAffineCaster->SetInput(  m_RigidOutputImage );
	movingAffineCaster->SetInput( downsampledAtlas );

  affineRegistration->SetFixedImage(    fixedAffineCaster->GetOutput()    );
  affineRegistration->SetMovingImage(   movingAffineCaster->GetOutput()   );

  fixedAffineCaster->Update();

  affineRegistration->SetFixedImageRegion( 
       fixedAffineCaster->GetOutput()->GetBufferedRegion() );	

	// Start affine registration
	std::cout<<"Start Affine Registration..."<<std::endl;
	itk::TimeProbe clock2;
	clock2.Start();
  try
  {
    affineRegistration->Update();
    std::cout << "Optimizer stop condition: "
              << affineRegistration->GetOptimizer()->GetStopConditionDescription()
              << std::endl;
  }
  catch( itk::ExceptionObject & err )
  {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
  }
	clock2.Stop();
	std::cout << "Time: " << clock2.GetTotal() << std::endl;

	AffineOptimizerType::ParametersType affineParameters = 
                  affineRegistration->GetLastTransformParameters();


	m_AffineTransform = AffineTransformType::New();
	m_AffineTransform->SetParameters(affineParameters);
	m_AffineTransform->SetCenter(affineTransform->GetCenter());

	// Resample the atlas image using the resulting affine transform
	typename ResampleFilterType::Pointer affineResampler = ResampleFilterType::New();
  typename ResampleInterpolatorType::Pointer affineResampleInterpolator = ResampleInterpolatorType::New();

	affineResampler->SetTransform( m_AffineTransform );
	affineResampler->SetInput( m_Atlas );
	affineResampler->SetInterpolator( affineResampleInterpolator );

	affineResampler->SetSize(    m_Atlas->GetLargestPossibleRegion().GetSize() );
	affineResampler->SetOutputOrigin(  m_Atlas->GetOrigin() );
	affineResampler->SetOutputSpacing( m_Atlas->GetSpacing() );
	affineResampler->SetOutputDirection( m_Atlas->GetDirection() );
	affineResampler->SetDefaultPixelValue( 0 );
	affineResampler->Update();

	m_AffineOutputImage = affineResampler->GetOutput();

	// Resample the subject image using the resulting rigid transform using original resolution
	typename ResampleFilterType::Pointer hrRigidResampler = ResampleFilterType::New();

	hrRigidResampler->SetTransform( m_RigidTransform );
	hrRigidResampler->SetInput( m_InputImage );
	hrRigidResampler->SetInterpolator( rigidResampleInterpolator );

	hrRigidResampler->SetSize(    m_Atlas->GetLargestPossibleRegion().GetSize() );
	hrRigidResampler->SetOutputOrigin(  m_Atlas->GetOrigin() );
	hrRigidResampler->SetOutputSpacing( m_Atlas->GetSpacing() );
	hrRigidResampler->SetOutputDirection( m_Atlas->GetDirection() );
	hrRigidResampler->SetDefaultPixelValue( 0 );
	hrRigidResampler->Update();

	m_RigidOutputImage = hrRigidResampler->GetOutput();

	// Resample the atlas mask using the resulting affine transform
	typename ReaderType::Pointer maskReader = ReaderType::New();
	std::string maskName = m_AtlasPath + "/mni_icbm152_t1_tal_nlin_sym_09c_mask.nii";
	maskReader->SetFileName(maskName);
	try
	{
		maskReader->Update();
	}
	catch(itk::ExceptionObject & err)
	{
		std::cerr<< err <<std::endl;
	}


	typename ResampleFilterType::Pointer maskResampler = ResampleFilterType::New();
	maskResampler->SetTransform( m_AffineTransform );
	maskResampler->SetInput( maskReader->GetOutput() );
	maskResampler->SetInterpolator( affineResampleInterpolator );

	maskResampler->SetSize( m_Atlas->GetLargestPossibleRegion().GetSize() );
	maskResampler->SetOutputOrigin(  m_Atlas->GetOrigin() );
	maskResampler->SetOutputSpacing( m_Atlas->GetSpacing() );
	maskResampler->SetOutputDirection( m_Atlas->GetDirection() );
	maskResampler->SetDefaultPixelValue( 0 );
	maskResampler->Update();

	m_AffineMask = maskResampler->GetOutput();


}


template<class TImage>
void BrainAlignment<TImage>::UpdateRSS()
{
	const unsigned int Dimension = 3;
	typedef itk::RegularStepGradientDescentOptimizer       OptimizerType;

  typedef itk::MattesMutualInformationImageToImageMetric< 
                                    InternalImageType, 
                                    InternalImageType >    MetricType;
  
  typedef itk::LinearInterpolateImageFunction< 
                                    InternalImageType,
                                    double          >    InterpolatorType;

  typedef itk::MultiResolutionImageRegistrationMethod< 
                                    InternalImageType, 
                                    InternalImageType >   RegistrationType;

	typedef itk::ScaleSkewVersor3DTransform<double>					TransformType;

  typedef itk::MultiResolutionPyramidImageFilter<
                                    InternalImageType,
                                    InternalImageType >   FixedImagePyramidType;
  typedef itk::MultiResolutionPyramidImageFilter<
                                    InternalImageType,
                                    InternalImageType >   MovingImagePyramidType;

  typedef itk::CastImageFilter< 
                        ImageType, InternalImageType > FixedCastFilterType;
  typedef itk::CastImageFilter< 
                        ImageType, InternalImageType > MovingCastFilterType;

	typedef itk::ImageFileReader<ImageType>	ReaderType;

	typedef itk::ResampleImageFilter< 
                            ImageType, 
                            ImageType >    ResampleFilterType;

  typedef itk::LinearInterpolateImageFunction< 
                                 ImageType, double >  ResampleInterpolatorType;

	typedef itk::ImageFileReader<InternalImageType>	FloatReaderType;

	typedef itk::DualResampleImageFilter< 
                            InternalImageType, 
                            InternalImageType >    FloatResampleFilterType;

  typedef itk::LinearInterpolateImageFunction< 
                                 InternalImageType, double >  FloatResampleInterpolatorType;

	typedef itk::DualResampleImageFilter<ImageType, ImageType>	DualResampleFilterType;


	// Load the atlas image
	typename ReaderType::Pointer atlasReader = ReaderType::New();
	std::string atlasName = m_AtlasPath + "/mni_icbm152_t1_tal_nlin_sym_09c.nii";
	atlasReader->SetFileName(atlasName);
	try
	{
		atlasReader->Update();
	}
	catch(itk::ExceptionObject & err)
	{
		std::cerr<< err <<std::endl;
	}

	m_Atlas = atlasReader->GetOutput();

	// Downsample the atlas image
	typename ImageType::Pointer downsampledAtlas = DownsampleImage(m_Atlas);

  MetricType::Pointer         metric        = MetricType::New();
  OptimizerType::Pointer			optimizer     = OptimizerType::New();
  InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
  RegistrationType::Pointer   registration  = RegistrationType::New();

  FixedImagePyramidType::Pointer fixedRigidImagePyramid = 
      FixedImagePyramidType::New();
  MovingImagePyramidType::Pointer movingRigidImagePyramid =
      MovingImagePyramidType::New();

	// Key components of the rigid registration
  registration->SetMetric(        metric        );
  registration->SetOptimizer(     optimizer     );
  registration->SetInterpolator(  interpolator  );
  registration->SetFixedImagePyramid( fixedRigidImagePyramid );
  registration->SetMovingImagePyramid( movingRigidImagePyramid ); 

  TransformType::Pointer  transform = TransformType::New();
  registration->SetTransform( transform );

	// Set parameter for the MI metric
	const int numberOfPixels = m_Atlas->GetBufferedRegion().GetNumberOfPixels();
  metric->SetNumberOfHistogramBins(50);
  metric->SetNumberOfSpatialSamples(m_NumberOfSamples);

	// Load the brain mask
	typename ReaderType::Pointer maskReader = ReaderType::New();
	std::string maskName = m_AtlasPath + "/mni_icbm152_t1_tal_nlin_sym_09c_mask.nii";
	maskReader->SetFileName(maskName);
	try
	{
		maskReader->Update();
	}
	catch(itk::ExceptionObject & err)
	{
		std::cerr<< err <<std::endl;
	}

	// Cast the short image type to float image type
  typename FixedCastFilterType::Pointer fixedCaster   = FixedCastFilterType::New();
  typename MovingCastFilterType::Pointer movingCaster = MovingCastFilterType::New();

	movingCaster->SetInput(  downsampledAtlas );
	fixedCaster->SetInput( m_InputImage );

  registration->SetFixedImage(    fixedCaster->GetOutput()    );
  registration->SetMovingImage(   movingCaster->GetOutput()   );

  fixedCaster->Update();

  registration->SetFixedImageRegion( 
       fixedCaster->GetOutput()->GetBufferedRegion() );	

  registration->SetInitialTransformParameters( transform->GetParameters() );

  typedef OptimizerType::ScalesType		OptimizerScalesType;
  OptimizerScalesType optimizerScales(	transform->GetNumberOfParameters() );
  optimizerScales[0] =  1.0; //Versor X scale
  optimizerScales[1] =  1.0; //Versor Y scale
  optimizerScales[2] =  1.0; //Versor Z scale
  optimizerScales[3] =  0.001; //Translation X scale
  optimizerScales[4] =  0.001; //Translation Y scale
  optimizerScales[5] =  0.001; //Translation Z scale
  optimizerScales[6] =  10.0;  //Scale X scale
  optimizerScales[7] =  10.0;  //Scale Y scale
  optimizerScales[8] =  10.0;  //Scale Z scale
  optimizerScales[9] =  10.0;  // [9] to [14] are 6 skew scales
  optimizerScales[10] =  1000.0;
  optimizerScales[11] =  1000.0;
  optimizerScales[12] =  1000.0;
  optimizerScales[13] =  1000.0;
  optimizerScales[14] =  1000.0;
  optimizer->SetScales( optimizerScales );

  optimizer->SetMaximumStepLength( 1 );
  optimizer->SetMinimumStepLength( 0.01 );
  optimizer->SetNumberOfIterations( 1200 );

	optimizer->MinimizeOn();

	// Start RSS registration
	std::cout<<"Start RigidScaleSkew Registration..."<<std::endl;
	itk::TimeProbe clock2;
	clock2.Start();
  try
  {
    registration->Update();
    std::cout << "Optimizer stop condition: "
              << registration->GetOptimizer()->GetStopConditionDescription()
              << std::endl;
  }
  catch( itk::ExceptionObject & err )
  {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
  }
	clock2.Stop();
	std::cout << "Time: " << clock2.GetTotal() << std::endl;

	OptimizerType::ParametersType parameters = 
                  registration->GetLastTransformParameters();


	m_RSSTransform = TransformType::New();
	m_RSSTransform->SetParameters(parameters);
	m_RSSTransform->SetCenter(transform->GetCenter());

	// Resample the atlas image using the resulting affine transform
	typename ResampleFilterType::Pointer rssResampler = ResampleFilterType::New();
    typename ResampleInterpolatorType::Pointer rssResampleInterpolator = ResampleInterpolatorType::New();

	rssResampler->SetTransform( m_RSSTransform );
	rssResampler->SetInput( maskReader->GetOutput() );
	rssResampler->SetInterpolator( rssResampleInterpolator );

	rssResampler->SetSize(    m_InputImage->GetLargestPossibleRegion().GetSize() );
	rssResampler->SetOutputOrigin(  m_InputImage->GetOrigin() );
	rssResampler->SetOutputSpacing( m_InputImage->GetSpacing() );
	rssResampler->SetOutputDirection( m_InputImage->GetDirection() );
	rssResampler->SetDefaultPixelValue( 0 );
	rssResampler->Update();

	m_RSSOutputImage = rssResampler->GetOutput();
}

#endif
