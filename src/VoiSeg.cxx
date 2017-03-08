#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkImageRegionIterator.h"

#include "itkVector.h"
#include "itkListSample.h"

#include "itkGaussianMixtureModelComponent.h"
#include "itkExpectationMaximizationMixtureModelEstimator.h"

#include "itkNormalVariateGenerator.h"

#include <itkMRIBiasFieldCorrectionFilter.h>

#include "itkThresholdImageFilter.h"

#include "itkMinimumMaximumImageCalculator.h"

#include "itkN4MRIBiasFieldCorrectionImageFilter.h"

typedef itk::Array< double > ParametersType;

const double PI = 4.0*atan(1.0);

template <class TImage>
void  EMSegmentation(TImage* image,const unsigned int numberOfClasses, std::vector< ParametersType > &initialParameters,itk::Array< double > &initialProportions, std::vector<ParametersType> &finalParameters, itk::Array< double > &finalProportion)
{
	const unsigned int vectorSize = 1;
	
	typedef itk::Vector< double, vectorSize > MeasurementVectorType;
	typedef itk::Statistics::ListSample< MeasurementVectorType > SampleType;
	SampleType::Pointer sample = SampleType::New();
	sample->SetMeasurementVectorSize( vectorSize ); // length of measurement vectors

	// Put image pixels into sample
	typedef TImage::PixelType PixelType;
	PixelType* image_pointer = image->GetPixelContainer()->GetImportPointer();

	const unsigned int image_dimension = TImage::ImageDimension;
	unsigned int image_length = 1;
	for(int i=0;i<image_dimension;i++)
	{
		image_length *= image->GetBufferedRegion().GetSize()[i];
	}

	for(int i=0;i<image_length;i++)
	{
		//if( image_pointer[i] != 0 )
			sample->PushBack(image_pointer[i]);
	}


	typedef itk::Statistics::GaussianMixtureModelComponent< SampleType > ComponentType;

	std::vector< ComponentType::Pointer > components;
	for ( unsigned int i = 0 ; i < numberOfClasses ; i++ )
	{
		components.push_back( ComponentType::New() );
		(components[i])->SetSample( sample );
		(components[i])->SetParameters( initialParameters[i] );
	}

	typedef itk::Statistics::ExpectationMaximizationMixtureModelEstimator< 
						   SampleType > EstimatorType;
	EstimatorType::Pointer estimator = EstimatorType::New();

	estimator->SetSample( sample );
	estimator->SetMaximumIteration( 200 );


	estimator->SetInitialProportions( initialProportions );

	for ( unsigned int i = 0 ; i < numberOfClasses ; i++)
	{
		estimator->AddComponent( (ComponentType::Superclass*)
								 (components[i]).GetPointer() );
	}

	estimator->Update();

	for ( unsigned int i = 0 ; i < numberOfClasses ; i++ )
	{
		std::cout << "Cluster[" << i << "]" << std::endl;
		std::cout << "    Parameters:" << std::endl;
		std::cout << "         " << (components[i])->GetFullParameters() 
				  << std::endl;
		std::cout << "    Proportion: ";
		std::cout << "         " << (estimator->GetProportions())[i] << std::endl;

		finalParameters.push_back((components[i])->GetFullParameters());
		
	}

	finalProportion = (estimator->GetProportions());

}

double ClassProbability(short p, ParametersType para, double proportion)
{
	float mean = para[0];
	float var = (para[1]);

	float d = p - mean;

	double amp = proportion/sqrt(2*PI*var);

	return amp*exp(-0.5*(d*d)/var);
}

int main(int argc, char* argv[])
{
 	if(argc < 4)
	{
		std::cout<<"Usage: "<<argv[0]<<" inputImage inputParameter GMMParameter"<<std::endl;
		return EXIT_FAILURE;
	}

	// Read input image
	typedef signed short PixelType;
	const unsigned int ImageDimension = 3;
	typedef itk::Image<PixelType, ImageDimension>	ImageType; 
	typedef itk::ImageFileReader<ImageType>	ReaderType;
	typedef itk::ImageFileWriter<ImageType> WriterType;

	std::cout<<"Start reading image: "<<std::flush;

	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(argv[1]);
	try
	{
		reader->Update();
	}
	catch(itk::ExceptionObject &err)
	{
		std::cerr<<err<<std::endl;
		return EXIT_FAILURE;
	}
	std::cout<<" [OK]"<<std::endl;

	ImageType::Pointer inputImage = reader->GetOutput();

	ImageType::SpacingType inputSpacing = inputImage->GetSpacing();

	// Read parameters
	std::ifstream infile;
	infile.open(argv[2]);
	if(infile.fail())
	{
		std::cerr<<"Cannot find parameter file."<<std::endl;
		return EXIT_FAILURE;
	}
	
	int directionIndex[3];

	ImageType::DirectionType inputDirection = inputImage->GetDirection();
	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			if(inputDirection[i][j] == 1 || inputDirection[i][j] == -1)
				directionIndex[i] = j;
		}
	}
	
	ImageType::PointType voiOrigin;
	ImageType::SizeType voiSize;
	for(int i=0;i<3;i++)
	{
		infile >> 	voiOrigin[i];
	}
	infile >> voiSize[2];
	infile >> voiSize[1];
	infile >> voiSize[0];
	
	//voiSize[0] = 20;
	//voiSize[1] = 20;
	//voiSize[2] = 20;

	// Create voi image
	std::cout<<"Start reading image: "<<std::flush;

	//ImageType::PointType voiOrigin;
	//voiOrigin[0] = 13.5;
	//voiOrigin[1] = 60.6;
	//voiOrigin[2] = -22.5;
	//ImageType::SizeType voiSize;
	//voiSize[0] = 20;
	//voiSize[1] = 15;
	//voiSize[2] = 20;

	ImageType::IndexType voiCenterIndex;
	inputImage->TransformPhysicalPointToIndex(voiOrigin, voiCenterIndex);

	ImageType::PixelType p = inputImage->GetPixel(voiCenterIndex);
	
	ImageType::PointType voiOriginPoint;
	//voiOriginPoint[0] = voiOrigin[0] + voiSize[0]/2;
	//voiOriginPoint[1] = voiOrigin[1] - voiSize[1]/2;
	//voiOriginPoint[2] = voiOrigin[2] - voiSize[2]/2;
	ImageType::IndexType voiOriginIndex;
	//inputImage->TransformPhysicalPointToIndex(voiOriginPoint, voiOriginIndex);
	voiOriginIndex[0] = voiCenterIndex[0] - voiSize[0]/2;
	voiOriginIndex[1] = voiCenterIndex[1] - voiSize[1]/2;
	voiOriginIndex[2] = voiCenterIndex[2] - voiSize[2]/2;

	ImageType::Pointer voiImage = ImageType::New();
	voiImage->CopyInformation(inputImage);
	ImageType::RegionType voiRegion;
	voiRegion.SetIndex(voiOriginIndex);
	voiRegion.SetSize(voiSize);
	voiImage->SetRegions(voiRegion);
	voiImage->Allocate();

	typedef itk::ImageRegionIterator<ImageType>	IteratorType;
	IteratorType inputIt(inputImage, voiRegion);
	IteratorType voiIt(voiImage, voiRegion);

	for(inputIt.GoToBegin(),voiIt.GoToBegin();!inputIt.IsAtEnd();++inputIt,++voiIt)
	{
		voiIt.Set(inputIt.Get());
	}

	WriterType::Pointer voiWriter = WriterType::New();
	voiWriter->SetInput(voiImage);
	voiWriter->SetFileName("voi.nii");
	voiWriter->Update();


	// Rescale the image to 0 - 255
	std::cout<<"Start intensity rescaling: "<<std::flush;
	typedef itk::RescaleIntensityImageFilter<ImageType, ImageType>	RescalerType;
	RescalerType::Pointer rescaler = RescalerType::New();
	rescaler->SetInput(voiImage);
	rescaler->SetOutputMaximum(255);
	rescaler->SetOutputMinimum(0);
	rescaler->Update();
	std::cout<<" [OK]"<<std::endl;

	ImageType::Pointer image = rescaler->GetOutput();


	// Bias Fied Correction
	std::cout<<"Start bias field correction: "<<std::flush;
	typedef itk::N4MRIBiasFieldCorrectionImageFilter<ImageType, ImageType> CorrecterType;
	CorrecterType::Pointer correcter=CorrecterType::New();
	correcter->SetInput(image);
	correcter->Update();
	std::cout<<" [OK]"<<std::endl;

	ImageType::Pointer correctedImage = correcter->GetOutput();

	WriterType::Pointer tempWriter = WriterType::New();
	tempWriter->SetInput(correctedImage);
	tempWriter->SetFileName("rescaledVoi.nii");
	tempWriter->Update();


	// Initial parameters read from txt file
	std::ifstream if_para(argv[3]);
	if(if_para.fail())
	{
		std::cerr<<"Cannot find parameter file."<<std::endl;
		return EXIT_FAILURE;
	}

	const unsigned int numberOfClasses = 3;

	std::vector< ParametersType > initialParameters( numberOfClasses );
	itk::Array< double > initialProportions(numberOfClasses);

	ParametersType params( numberOfClasses );
	for(int i=0;i<numberOfClasses;i++)
	{
		if_para >> 	params[0] >> params[1] >> initialProportions[i];
		initialParameters[i] = params;
	}


	std::vector<ParametersType> finalParameters;
	itk::Array< double > finalProportion;

	std::cout<<"Start EM estimation: "<<std::flush;
	EMSegmentation<ImageType>(image, numberOfClasses, initialParameters, initialProportions, finalParameters, finalProportion);
	std::cout<<" [OK]"<<std::endl;

	// Create labelled image using the final parameters
	ImageType::Pointer labelledImage = ImageType::New();
	labelledImage->CopyInformation(image);
	labelledImage->SetRegions(image->GetBufferedRegion());
	labelledImage->Allocate();
	labelledImage->FillBuffer(-100);

	int count[3] = {0,0,0};
	double pCount[3] = {0.0,0.0,0.0};

	typedef itk::ImageRegionIterator<ImageType> IteratorType;
	IteratorType imageIt(image, image->GetBufferedRegion());
	IteratorType labelledIt(labelledImage, labelledImage->GetBufferedRegion());
	for(imageIt.GoToBegin(),labelledIt.GoToBegin();!imageIt.IsAtEnd();++imageIt,++labelledIt)
	{
		short p = imageIt.Get();

		double c[numberOfClasses];
		double pt = 0.0;
		for(int i=0;i<numberOfClasses;i++)
		{
			c[i] = ClassProbability(p, finalParameters[i], finalProportion[i]);
			pt += c[i];
		}
		for(int i=0;i<numberOfClasses;i++)
		{
			pCount[i] += c[i] / pt;
		}

		int maxIdx = 0;
		double max = 0.0;
		for(int i=0;i<numberOfClasses;i++)
		{
			if(c[i]>max)
			{
				maxIdx = i;
				max = c[i];
			}
		}

		count[maxIdx]++;
		labelledIt.Set(maxIdx*100);
	}

	WriterType::Pointer writer = WriterType::New();
	writer->SetInput(labelledImage);
	writer->SetFileName("SegmentedVoi.nii");
	writer->Update();

	std::cout<<std::endl;
	std::cout<<"Results by point counting: "<<std::endl;
	int totalCount = voiSize[0] * voiSize[1] * voiSize[2];
	std::cout<<"WM: "<<(float)count[0]/totalCount<<std::endl;
	std::cout<<"GM: "<<(float)count[1]/totalCount<<std::endl;
	std::cout<<"CSF:"<<(float)count[2]/totalCount<<std::endl;
	std::cout<<std::endl;

	std::cout<<"Results by partial volume probability counting: "<<std::endl;
	float totalPCount = pCount[0] + pCount[1] + pCount[2];
	std::cout<<"WM: "<<(float)pCount[0]/totalCount<<std::endl;
	std::cout<<"GM: "<<(float)pCount[1]/totalCount<<std::endl;
	std::cout<<"CSF:"<<(float)pCount[2]/totalCount<<std::endl;


	return EXIT_SUCCESS;
}