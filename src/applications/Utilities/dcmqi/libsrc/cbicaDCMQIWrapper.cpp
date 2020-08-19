#include "dcmqi/Helper.h"
#include "dcmqi/Exceptions.h"
#include "dcmqi/ImageSEGConverter.h"
#include "json.h"
#include "itkImage.h"
#include "itkImageFileReader.h"

namespace cbica
{
  void ConvertNiftiToDicomSeg(const std::string& niftiSegmentationFile, const std::string& dicomReferenceDirectory, const std::string& dicomSegMetaJSON, const std::string& outputDicomFile)
  {
    using ShortImageType = itk::Image< short, 3 >; // dicom seg only takes short, apparently
    using ShortImageReaderType = itk::ImageFileReader< ShortImageType >;
    auto reader = ShortImageReaderType::New();
    reader->SetFileName(niftiSegmentationFile);
    reader->Update();
    auto outputSEGFileName = outputDicomFile;

    using helper = dcmqi::Helper;

    auto referenceDicomFiles = helper::getFileListRecursively(dicomReferenceDirectory.c_str());

    auto dcmDatasets = helper::loadDatasets(referenceDicomFiles);

    if (dcmDatasets.empty())
    {
      cerr << "Error: no DICOM could be loaded from the specified list/directory.\n";
      exit(EXIT_FAILURE);
    }

    std::ifstream metainfoStream(dicomSegMetaJSON.c_str(), std::ios_base::binary);
    std::string metadata((std::istreambuf_iterator<char>(metainfoStream)),
      (std::istreambuf_iterator<char>()));

    Json::Value metaRoot;
    istringstream metainfoisstream(metadata);
    metainfoisstream >> metaRoot;

    std::vector< typename ShortImageType::Pointer > segmentations; // we will only have one segmentation file
    segmentations.push_back(reader->GetOutput());

    if (metaRoot.isMember("segmentAttributesFileMapping"))
    {
      if (metaRoot["segmentAttributesFileMapping"].size() != metaRoot["segmentAttributes"].size())
      {
        cerr << "Number of files in segmentAttributesFileMapping should match the number of entries in segmentAttributes!\n";
        exit(EXIT_FAILURE);
      }

      // we only accept a single nifti segmentation object
      auto segImageFiles = std::vector< std::string >{ niftiSegmentationFile }; // this to have consistency with dcmqi api

      // otherwise, re-order the entries in the segmentAtrributes list to match the order of files in segmentAttributesFileMapping
      Json::Value reorderedSegmentAttributes;
      std::vector<int> fileOrder(segImageFiles.size());
      fill(fileOrder.begin(), fileOrder.end(), -1);
      std::vector< typename ShortImageType::Pointer > segmentationsReordered(segImageFiles.size());
      for (size_t filePosition = 0; filePosition < segImageFiles.size(); filePosition++)
      {
        for (size_t mappingPosition = 0; mappingPosition < segImageFiles.size(); mappingPosition++)
        {
          string mappingItem = metaRoot["segmentAttributesFileMapping"][static_cast<int>(mappingPosition)].asCString();
          size_t foundPos = segImageFiles[filePosition].rfind(mappingItem);
          if (foundPos != std::string::npos)
          {
            fileOrder[filePosition] = mappingPosition;
            break;
          }
        }
        if (fileOrder[filePosition] == -1)
        {
          cerr << "Failed to map " << segImageFiles[filePosition] << " from the segmentAttributesFileMapping attribute to an input file name!" << endl;
          exit(EXIT_FAILURE);
        }
      }
      cout << "Order of input ITK images updated as shown below based on the segmentAttributesFileMapping attribute:" << endl;
      for (size_t i = 0; i < segImageFiles.size(); i++)
      {
        cout << " image " << i << " moved to position " << fileOrder[i] << "\n";
        segmentationsReordered[fileOrder[i]] = segmentations[i];
      }
      segmentations = segmentationsReordered;
    }

    try
    {
      bool skipEmptySlices = false;
      DcmDataset* result = dcmqi::ImageSEGConverter::itkimage2dcmSegmentation(dcmDatasets, segmentations, metadata, skipEmptySlices);

      if (result == NULL)
      {
        std::cerr << "ERROR: Conversion failed." << std::endl;
        exit(EXIT_FAILURE);
      }
      else
      {
        DcmFileFormat segdocFF(result);
        bool compress = false;
        if (compress)
        {
          CHECK_COND(segdocFF.saveFile(outputSEGFileName.c_str(), EXS_DeflatedLittleEndianExplicit));
        }
        else
        {
          CHECK_COND(segdocFF.saveFile(outputSEGFileName.c_str(), EXS_LittleEndianExplicit));
        }

        std::cout << "Saved segmentation as " << outputSEGFileName << endl;
      }

      for (size_t i = 0; i < dcmDatasets.size(); i++)
      {
        delete dcmDatasets[i];
      }
      if (result != NULL)
        delete result;
    }
    catch (std::exception e)
    {
      std::cerr << "Fatal error encountered: " << e.what() << "\n";
    }
  }

};