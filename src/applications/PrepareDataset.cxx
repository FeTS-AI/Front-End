#include <array>
#include <stdio.h>
#include <iostream>
#include <fstream>

#include "cbicaCmdParser.h"
#include "cbicaUtilities.h"
#include "CaPTkGUIUtils.h"
#include "cbicaITKImageInfo.h"
#include "cbicaITKSafeImageIO.h"

#include "BiasCorrection.hpp"

std::vector< std::map< std::string, std::string > > GetCSVContents(const std::string &fileName)
{
  std::vector< std::map< std::string, std::string > > csvContents;
  std::ifstream data(fileName.c_str());
  std::string line, cell;
  std::vector< std::string > headers; // csv headers

  size_t i = 0, j = 0;
  while (std::getline(data, line))
  {
    j = 0;
    std::map< std::string, std::string > currentRow;
    std::stringstream lineStream(line);
    while (std::getline(lineStream, cell, ','))
    {
      if (i == 0)
      {
        auto temp = cell;
        std::transform(temp.begin(), temp.end(), temp.begin(), ::tolower);
        temp = cbica::stringReplace(temp, "-", ""); // remove all hyphens
        temp = cbica::stringReplace(temp, "_", ""); // remove all underscores

        if ((temp == "patientid") || (temp == "subjectid") || (temp == "subject") || (temp == "subid"))
        {
          headers.push_back("ID");
        }
        else if ((temp == "t1gd") || (temp == "t1ce") || (temp == "t1post"))
        {
          headers.push_back("T1GD");
        }
        else if ((temp == "t1") || (temp == "t1pre"))
        {
          headers.push_back("T1");
        }
        else if (temp == "t2")
        {
          headers.push_back("T2");
        }
        else if ((temp == "t2flair") || (temp == "flair") || (temp == "fl") || (temp.find("fl") != std::string::npos) || (temp.find("t2fl") != std::string::npos))
        { 
          headers.push_back("FLAIR");
        }
      }
      else
      {
        if (headers.size() != 5)
        {
          std::cerr << "All required headers were not found in CSV. Please ensure the following are present: 'PatientID,T1,T1GD,T2,T2FLAIR'";
          return csvContents;
        }
        if (cell.find(" ") != std::string::npos)
        {
          std::cerr << "Please ensure that there are no spaces in the file paths.";
          return csvContents;
        }
        //auto temp = cbica::stringReplace(cell, "\"", ""); // remove all double quotes
        //temp = cbica::stringReplace(temp, "'", ""); // remove all single quotes
        currentRow[headers[j]] = cell;
      }
      j++;
    }
    if (i != 0)
    {
      csvContents.push_back(currentRow);
    }
    i++;
  }

  return csvContents;
}

std::string getStdoutFromCommand(const std::string command)
{
  std::array<char, 256> buffer;
  std::string result;
#ifdef WIN32
#define PCLOSE _pclose
#define POPEN _popen
#else
#define PCLOSE pclose
#define POPEN popen
#endif
  std::unique_ptr<FILE, decltype(&PCLOSE)> pipe(POPEN(command.c_str(), "r"), PCLOSE);
  if (!pipe)
  {
    throw std::runtime_error("popen() failed!");
  }
  while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr)
  {
    result += buffer.data();
  }
  return result;
}

bool BraTSPipeline(const std::map< std::string, std::string >& inputFiles, const std::string& outputDir)
{
  auto debug = false;
  // sanity checks
  for (auto it = inputFiles.begin(); it != inputFiles.end(); it++)
  {
    if (!cbica::exists(it->second))
    {
      std::cerr << "Couldn't find the modality '" << it->first << "', denoted by '" << it->second << "'.\n";
      return false;
    }

    auto inputImageInfo = cbica::ImageInfo(it->second);
    if (inputImageInfo.GetImageDimensions() != 3)
    {
      std::cerr << "The BraTS pipeline is only valid for 3D images, whereas the image '"
        << it->second << "' for modality '" << it->first << "' has " <<
        inputImageInfo.GetImageDimensions() << " dimentions.\n";
      return false;
    }
  }
  using ImageType = itk::Image< float, 3 >; // default image type

  // variables to store various images
  std::map< std::string, ImageType::Pointer > inputImages, inputImages_processed;
  std::map< std::string, std::string > inputModalities_orientation;

  // default names
  std::map< std::string, std::string > outputNames,
    inputReorientedFiles, inputReorientedBiasFiles, // filenames for reoriented and bias-corrected files
    outputMatFiles, outputRegisteredImages, outputRegisteredMaskedImages; // filenames for  matrices and images

  if (debug)
  {
    std::cout << "Reading input images.\n";
  }

  // construct path to dcm2niix for debug/release modes and different OS
  std::string m_exe;
#ifdef CAPTK_PACKAGE_PROJECT
#if WIN32
  m_exe = cbica::getExecutablePath() + "/dcm2niix.exe";
#else
  m_exe = cbica::getExecutablePath() + "/dcm2niix";
#endif
#else
#if WIN32
  m_exe = std::string(PROJECT_SOURCE_DIR) + "/src/applications/individualApps/dcm2niix/dcm2niix.exe";
#else
  m_exe = std::string(PROJECT_SOURCE_DIR) + "/src/applications/individualApps/dcm2niix/dcm2niix";
#endif
#endif

  for (auto it = inputFiles.begin(); it != inputFiles.end(); it++)
  {
    auto modality = it->first;
    /// [1] read image - DICOM to NIfTI conversion, if applicable
    inputImages[modality] = cbica::ReadImage< ImageType >(it->second);

    if (inputImages[modality].IsNull() && cbica::IsDicom(it->second))
    {
      auto dicomFolderPath = cbica::getFilenamePath(it->second);

      if (!cbica::isFile(m_exe))
      {
        std::cerr << "Couldn't find the dcm2niix executable, which was expected in '" << cbica::normalizePath(m_exe) << "'.\n";
        return false;
      }
      //else
      //{

      //}

      auto tempOutput = cbica::createTemporaryDirectory();
      //construct command
      std::string fullCommandToRun = cbica::normPath(m_exe) + " -o " + cbica::normPath(tempOutput) + " -z y \"" + cbica::normPath(dicomFolderPath) + "\"";

      //run command via system call
      if (std::system((fullCommandToRun).c_str()) != 0)
      {
        std::cerr << "Something went wrong during dicom to nifti conversion, please re-try or contact sofware@cbica.upenn.edu.\n";
        return false;
      }

      auto filesInDir = cbica::filesInDirectory(tempOutput);
      for (size_t i = 0; i < filesInDir.size(); i++)
      {
        if (cbica::getFilenameExtension(filesInDir[i]) == ".nii.gz")
        {
          inputImages[modality] = cbica::ReadImage< ImageType >(filesInDir[i]);
          break;
        }
      }
      cbica::removeDirectoryRecursively(tempOutput, true);
    }

    if (inputImages[modality].IsNotNull())
    {
      auto fileToWrite = outputDir + "/raw_" + modality + ".nii.gz";
      cbica::WriteImage< ImageType >(inputImages[modality], fileToWrite);
    }
    else
    {
      if (cbica::IsDicom(it->second))
      {
        std::cerr << "Something went wrong with the DICOM to NIfTI conversion for modality '" <<
          modality << "' with filename '" << it->second << "'"
          << ", please use another package to conver to NIfTI and try again.\n";
        return false;
      }
      else
      {
        std::cerr << "Something went wrong with reading the raw input image, please re-try or contact sofware@cbica.upenn.edu.\n";
        return false;
      }
    }

    /// [2] LPS/RAI re-orientation
    if (debug)
    {
      std::cout << "Performing re-orientation to LPS for modality '" << modality << "'.\n";
    }

    auto temp = cbica::GetImageOrientation(inputImages[modality], "RAI");
    inputModalities_orientation[modality] = temp.first;
    inputImages_processed[modality] = temp.second;
    if (inputImages_processed[modality].IsNull())
    {
      std::cerr << "Something went wrong with re-orienting the input image, please re-try or contact sofware@cbica.upenn.edu.\n";
      return false;
    }
    else
    {
      inputReorientedFiles[modality] = outputDir + "/raw_rai_" + modality + ".nii.gz";
      // the re-oriented images need to be written because these are passed on to greedy
      cbica::WriteImage< ImageType >(inputImages_processed[modality], inputReorientedFiles[modality]);
    }

    /// [3] N4 bias correction

    if (debug)
    {
      std::cout << "Starting bias correction for modality '" << modality << "'.\n";
    }

    // the bias-corrected images need to be written because these are passed on to greedy
    inputReorientedBiasFiles[modality] = outputDir + "/raw_rai_n4_" + modality + ".nii.gz";

    if (!cbica::fileExists(inputReorientedBiasFiles[modality]))
    {
      BiasCorrection biasCorrector;
      {
        using MaskImageType = itk::Image<unsigned char, ImageType::ImageDimension>;
        typename MaskImageType::Pointer maskImage; // mask inits to null
        auto outputImage = biasCorrector.Run<TImageType, MaskImageType>("n4",
          inputImages_processed[modality],
          maskImage,
          BiasCorrection::default_splineOrder,
          BiasCorrection::default_maxIterations,
          BiasCorrection::default_fittingLevels,
          BiasCorrection::default_filterNoise,
          BiasCorrection::default_fwhm,
          BiasCorrection::default_otsuBins);
        if (outputImage.IsNotNull())
        {
          inputImages_processed[modality] = outputImage;
          inputImages_processed[modality]->DisconnectPipeline();
        }
        else
        {
          std::cerr << "Something went wrong with bias-correcting the re-oriented image, please re-try or contact sofware@cbica.upenn.edu.\n";
          return false;
        }
      }
      cbica::WriteImage< ImageType >(inputImages_processed[modality], inputReorientedBiasFiles[modality]);
    }

    if (modality != "T1CE")
    {
      outputNames[modality] = modality + "_to_T1CE"; // all output names can be controlled from here
    }
    else
    {
      outputNames[modality] = modality + "_to_SRI"; // all output names can be controlled from here
    }
  } // end inputFiles iterator

  /// [4] Registration using Greedy
  // we do T1CE to Atlas registration first because other registrations are dependent on this
  if (debug)
  {
    std::cout << "Registering T1CE to SRI atlas.\n";
  }

  auto greedyPathAndDim = getApplicationPath("greedy") + " -d 3";

  auto captkDataDir = getCaPTkDataDir();
  auto atlasImage = captkDataDir + "/sri24/atlasImage.nii.gz";
  outputMatFiles["T1CE"] = outputDir + "/" + outputNames["T1CE"] + ".mat";
  outputRegisteredImages["T1CE"] = outputDir + "/" + outputNames["T1CE"] + ".nii.gz";
  outputRegisteredMaskedImages["T1CE"] = outputDir + "/brain_T1CE.nii.gz";

  std::string fullCommand;

  if (!cbica::exists(outputMatFiles["T1CE"]))
  {
    fullCommand = " -a -m NMI -i " + atlasImage + " " + inputReorientedBiasFiles["T1CE"]
      + " -o " + outputMatFiles["T1CE"] + " -ia-image-centers -n 100x50x10 -dof 6";

    if (debug)
    {
      std::cout << "Greedy command: " << greedyPathAndDim + fullCommand << "\n";
    }
    if (std::system((greedyPathAndDim + fullCommand).c_str()) != 0)
    {
      std::cerr << "Something went wrong when registering T1CE image to SRI atlas, please re-try or contact sofware@cbica.upenn.edu.\n";
      return false;
    }
  } // end outputMatFiles["T1CE"] check

  if (!cbica::exists(outputRegisteredImages["T1CE"]))
  {
    fullCommand = " -rf " + atlasImage + " -ri LINEAR -rm " +
      inputReorientedFiles["T1CE"] + " " + outputRegisteredImages["T1CE"] + " -r " +
      outputMatFiles["T1CE"];

    if (debug)
    {
      std::cout << "Greedy command: " << greedyPathAndDim + fullCommand << "\n";
    }

    if (std::system((greedyPathAndDim + fullCommand).c_str()) != 0)
    {
      std::cerr << "Something went wrong when applying registration matrix to generate T1CE image in SRI atlas space, please re-try or contact sofware@cbica.upenn.edu.\n";
      return false;
    }
  } // end outputRegisteredImages["T1CE"] check

  for (auto it = inputFiles.begin(); it != inputFiles.end(); it++)
  {
    auto modality = it->first;
    if (modality != "T1CE") // T1CE registration has happened before
    {
      outputMatFiles[modality] = outputDir + "/" + outputNames[modality] + ".mat";
      outputRegisteredImages[modality] = outputDir + "/" + modality + "_to_SRI.nii.gz";
      outputRegisteredMaskedImages[modality] = outputDir + "/brain_" + modality + ".nii.gz";

      if (!cbica::exists(outputMatFiles[modality]))
      {
        // we use the bias-corrected image for registration as it is easier localize transformations
        fullCommand = " -a -m NMI -i " + inputReorientedBiasFiles["T1CE"] + " " + inputReorientedBiasFiles[modality]
          + " -o " + outputMatFiles[modality] + " -ia-image-centers -n 100x50x10 -dof 6";
        if (debug)
        {
          std::cout << "Registering " << modality << " to T1CE.\n";
          std::cout << "Greedy command: " << greedyPathAndDim + fullCommand << "\n";
        }

        if (std::system((greedyPathAndDim + fullCommand).c_str()) != 0)
        {
          std::cerr << "Something went wrong when registering " << modality
            << "to T1CE image, please re-try or contact sofware@cbica.upenn.edu.\n";
          return false;
        }
      } // end outputMatFiles[modality] check

      if (debug)
      {
        std::cout << "Generating image for " << modality << " registered to the atlas.\n";
        std::cout << "Greedy command: " << greedyPathAndDim + fullCommand << "\n";
      }

      if (!cbica::exists(outputRegisteredImages[modality]))
      {
        // the final registration is applied on the original image after re-orientation (not bias-corrected) to
        // ensure maximum fidelity with original image
        fullCommand = " -rf " + atlasImage + " -ri LINEAR -rm " + inputReorientedFiles[modality] + " " +
          outputRegisteredImages[modality] + " -r "
          + outputMatFiles["T1CE"] + " "
          + outputMatFiles[modality];

        if (std::system((greedyPathAndDim + fullCommand).c_str()) != 0)
        {
          std::cerr << "Something went wrong when applying registration matrix to generate " << modality << " image in SRI atlas space, please re-try or contact sofware@cbica.upenn.edu.\n";
          return false;
        }
      } // end outputRegisteredImages[modality] check
    } // end modality check
  } // end modality loop

  // variables that are used later on
  auto finalBrainMask = cbica::normalizePath(outputDir + "/brainMask_SRI.nii.gz");
  auto deepMedicExe = getApplicationPath("DeepMedic");
  auto brainMaskFile = outputDir + "/dmOut_skull/brainMask_SRI.nii.gz";

  /// [5] Skull-stripping using DeepMedic
  if (debug)
  {
    std::cout << "Starting skull-stripping using DeepMedic.\n";
  }

  if (!cbica::exists(brainMaskFile))
  {
    fullCommand = " -md " + captkDataDir + "/deepMedic/saved_models/skullStripping/ " +
      "-i " + outputRegisteredImages["T1"] + "," +
      outputRegisteredImages["T1CE"] + "," +
      outputRegisteredImages["T2"] + "," +
      outputRegisteredImages["FL"] + " -o " +
      brainMaskFile;

    if (debug)
    {
      std::cout << "Command for DeepMedic: " << deepMedicExe + fullCommand << "\n";
    }

    if (std::system((deepMedicExe + fullCommand).c_str()) != 0)
    {
      std::cerr << "Something went wrong when performing skull-stripping using DeepMedic, please re-try or contact sofware@cbica.upenn.edu.\n";
      return false;
    }
  } // end brainMask check

  if (!cbica::exists(brainMaskFile))
  {
    std::cerr << "Brain Mask was not written, cannot proceed.\n";
    return false;
  }

  // variables to store outputs in patient space
  std::map< std::string, std::string > outputFiles_withoutOrientationFix, outputFiles_withOrientationFix;

  cbica::WriteImage< TImageType >(
    cbica::ReadImage< TImageType >(brainMaskFile),
    finalBrainMask
    );

  // iterate over outputRegisteredMaskedImages
  for (auto it = outputRegisteredMaskedImages.begin(); it != outputRegisteredMaskedImages.end(); it++)
  {
    auto modality = it->first;
    auto maskFilter = itk::MaskImageFilter< ImageType, ImageType >::New();
    maskFilter->SetInput(cbica::ReadImage< ImageType >(outputRegisteredImages[modality]));
    maskFilter->SetMaskImage(cbica::ReadImage< TImageType >(finalBrainMask));
    try
    {
      maskFilter->Update();
    }
    catch (const std::exception& e)
    {
      std::cerr << "Something went wrong when applying the brain mask to modality '"
        << modality << "': " << e.what();
      return false;
    }
    cbica::WriteImage< ImageType >(maskFilter->GetOutput(), it->second); // write the masked image 
  }

}

int main(int argc, char** argv)
{
  cbica::CmdParser parser(argc, argv, "PrepareDataset");

  parser.addApplicationDescription("This application calls the BraTSPipeline for all input images and stores the final and intermediate files separately");
  parser.addRequiredParameter("i", "inputCSV", cbica::Parameter::FILE, "Input CSV file", "Input CSV file which contains paths to structural images", "Headers should be 'PatientID,T1,T1GD,T2,T2FLAIR'");
  parser.addRequiredParameter("o", "outputDir", cbica::Parameter::DIRECTORY, "Directory", "Output directory for final output", "This will write 2 folders: 'DataForFeTS' and 'DataForQC'", 
    "Former contains only the files needed for FeTS inference/training and ", "latter contains all intermediate files from this processing");

  std::string inputCSV, outputDir;

  parser.getParameterValue("i", inputCSV);
  parser.getParameterValue("o", outputDir);

  auto csvContents = GetCSVContents(inputCSV);

  if (csvContents.empty())
  {
    std::cerr << "Parsed CSV data structure is empty, cannot proceed.\n";
    return EXIT_FAILURE;
  }

  // set up the output directories
  auto outputDir_qc = cbica::normPath(outputDir + "/DataForQC");
  auto outputDir_final = cbica::normPath(outputDir + "/DataForFeTS");
  cbica::createDir(outputDir);
  cbica::createDir(outputDir_qc);
  cbica::createDir(outputDir_final);

  auto bratsPipeline_exe = getApplicationPath("BraTSPipeline");

  if (!cbica::isFile(bratsPipeline_exe))
  {
    std::cerr << "BraTSPipeline was not found in the installation, cannot proceed.\n";
    return EXIT_FAILURE;
  }
  // iterate through all subjects
  for (size_t i = 0; i < csvContents.size(); i++)
  {
    std::cout << "Started processing subject '" << csvContents[i]["ID"] << "'\n";

    auto interimOutputDir = outputDir_qc + "/" + csvContents[i]["ID"];
    auto finalSubjectOutputDir = outputDir_final + "/" + csvContents[i]["ID"];
    cbica::createDir(finalSubjectOutputDir);

    auto command = bratsPipeline_exe + " -t1 " + csvContents[i]["T1"] + " -t1c " + csvContents[i]["T1GD"] + " -t2 " + csvContents[i]["T2"] + " -fl " + csvContents[i]["FLAIR"] + " -o " + interimOutputDir + " -s 1";

    std::map< std::string, std::string > inputFiles;
    inputFiles["T1"] = csvContents[i]["T1"];
    inputFiles["T1CE"] = csvContents[i]["T1GD"];
    inputFiles["T2"] = csvContents[i]["T2"];
    inputFiles["FLAIR"] = csvContents[i]["FLAIR"];
    
    if (BraTSPipeline(inputFiles, interimOutputDir))
    {
      auto msg = "BraTSPipeline failed for subject " + csvContents[i]["ID"];
      if (cbica::isFile(interimOutputDir + "/brain_T1CE.nii.gz"))
      {
        cbica::copyFile(interimOutputDir + "/brain_T1CE.nii.gz", finalSubjectOutputDir + "/" + csvContents[i]["ID"] + "_brain_t1ce.nii.gz");
      }
      if (cbica::isFile(interimOutputDir + "/brain_T1GD.nii.gz"))
      {
        cbica::copyFile(interimOutputDir + "/brain_T1GD.nii.gz", finalSubjectOutputDir + "/" + csvContents[i]["ID"] + "_brain_t1ce.nii.gz");
      }
      if (cbica::isFile(interimOutputDir + "/brain_T1.nii.gz"))
      {
        cbica::copyFile(interimOutputDir + "/brain_T1.nii.gz", finalSubjectOutputDir + "/" + csvContents[i]["ID"] + "_brain_t1.nii.gz");
      }
      if (cbica::isFile(interimOutputDir + "/brain_T2.nii.gz"))
      {
        cbica::copyFile(interimOutputDir + "/brain_T2.nii.gz", finalSubjectOutputDir + "/" + csvContents[i]["ID"] + "_brain_t2.nii.gz");
      }
      if (cbica::isFile(interimOutputDir + "/brain_FL.nii.gz"))
      {
        cbica::copyFile(interimOutputDir + "/brain_FL.nii.gz", finalSubjectOutputDir + "/" + csvContents[i]["ID"] + "_brain_flair.nii.gz");
      }
    }

    //auto log = getStdoutFromCommand(command);
    //std::ofstream myfile;
    //myfile.open(interimOutputDir + "/log.txt");
    //myfile << log;
    //myfile.close();

    //auto msg = "BraTSPipeline failed for subject " + csvContents[i]["ID"];
    //if (cbica::isFile(interimOutputDir + "/brain_T1CE.nii.gz"))
    //{
    //  cbica::copyFile(interimOutputDir + "/brain_T1CE.nii.gz", finalSubjectOutputDir + "/" + csvContents[i]["ID"] + "_brain_t1ce.nii.gz");
    //}
    //else if (cbica::isFile(interimOutputDir + "/brain_T1GD.nii.gz"))
    //{
    //  cbica::copyFile(interimOutputDir + "/brain_T1GD.nii.gz", finalSubjectOutputDir + "/" + csvContents[i]["ID"] + "_brain_t1ce.nii.gz");
    //}
    //else
    //{
    //  std::cerr << msg << "\n";
    //}
    //if (cbica::isFile(interimOutputDir + "/brain_T1.nii.gz"))
    //{
    //  cbica::copyFile(interimOutputDir + "/brain_T1.nii.gz", finalSubjectOutputDir + "/" + csvContents[i]["ID"] + "_brain_t1.nii.gz");
    //}
    //else
    //{
    //  std::cerr << msg << "\n";
    //}
    //if (cbica::isFile(interimOutputDir + "/brain_T2.nii.gz"))
    //{
    //  cbica::copyFile(interimOutputDir + "/brain_T2.nii.gz", finalSubjectOutputDir + "/" + csvContents[i]["ID"] + "_brain_t2.nii.gz");
    //}
    //else
    //{
    //  std::cerr << msg << "\n";
    //}
    //if (cbica::isFile(interimOutputDir + "/brain_FL.nii.gz"))
    //{
    //  cbica::copyFile(interimOutputDir + "/brain_FL.nii.gz", finalSubjectOutputDir + "/" + csvContents[i]["ID"] + "_brain_flair.nii.gz");
    //}
    //else
    //{
    //  std::cerr << msg << "\n";
    //}
  }

  return EXIT_SUCCESS;
}