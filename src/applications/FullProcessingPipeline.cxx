#include <array>
#include <stdio.h>
#include <iostream>
#include <fstream>

#include "cbicaCmdParser.h"
#include "cbicaUtilities.h"
#include "CaPTkGUIUtils.h"
#include "cbicaITKImageInfo.h"
#include "cbicaITKSafeImageIO.h"


int main(int argc, char** argv)
{
  cbica::CmdParser parser(argc, argv, "FullProcessingPipeline");

  auto hardcodedNativeModelWeightPath = getCaPTkDataDir() + "/fets";
  auto allArchs = cbica::subdirectoriesInDirectory(hardcodedNativeModelWeightPath);
  std::string allArchsString;
  for (size_t i = 0; i < allArchs.size(); i++)
  {
    allArchsString += allArchs[i] + ",";
  }
  allArchsString.pop_back();

  std::string modelName, loggingDir, archs = "3dresunet", fusionMethod = "STAPLE";

  parser.addApplicationDescription("This application calls the BraTSPipeline for all input images and stores the final and intermediate files separately");
  parser.addRequiredParameter("i", "inputCSV", cbica::Parameter::FILE, "Input CSV file", "Input CSV file which contains paths to structural images", "Headers should be 'PatientID,T1,T1GD,T2,T2FLAIR'");
  parser.addRequiredParameter("o", "outputDir", cbica::Parameter::DIRECTORY, "Directory", "Output directory for final output", "This will write 2 folders: 'DataForFeTS' and 'DataForQC'", 
    "Former contains only the files needed for FeTS inference/training and ", "latter contains all intermediate files from this processing");
  parser.addOptionalParameter("a", "archs", cbica::Parameter::STRING, allArchsString, "The architecture(s) to infer/train on", "Only a single architecture is supported for training", "Comma-separated values for multiple options", "Defaults to: " + archs);
  parser.addOptionalParameter("lF", "labelFuse", cbica::Parameter::STRING, "STAPLE,ITKVoting,SIMPLE,MajorityVoting", "The label fusion strategy to follow for multi-arch inference", "Comma-separated values for multiple options", "Defaults to: " + fusionMethod);
  parser.addOptionalParameter("g", "gpu", cbica::Parameter::BOOLEAN, "0-1", "Whether to run the process on GPU or not", "Defaults to '0'");
  parser.addOptionalParameter("L", "LoggingDir", cbica::Parameter::DIRECTORY, "Dir with write access", "Location of logging directory");

  auto preparedataset_path = getApplicationPath("PrepareDataset");

  std::string inputCSV, outputDir, gpu_request = "0";

  // parser CLI parameters
  parser.getParameterValue("i", inputCSV);
  parser.getParameterValue("o", outputDir);

  if (parser.isPresent("L"))
  {
    parser.getParameterValue("L", loggingDir);
  }
  else
  {
    loggingDir = cbica::createTemporaryDirectory() + "/logs";
    std::cout << "Using the following directory as logging directory: " << loggingDir << "\n";
    cbica::createDir(loggingDir);
  }

  if (parser.isPresent("a"))
  {
    parser.getParameterValue("a", archs);
  }
  else
  {
    std::cerr << "Please specify at least 2 architectures on which to perform inference.\n";
  }
  if (parser.isPresent("lF"))
  {
    parser.getParameterValue("lF", fusionMethod);
  }
  if (parser.isPresent("g"))
  {
    parser.getParameterValue("g", gpu_request);
  }
  
  // run preparedataset
  auto full_command = preparedataset_path + " -i " + inputCSV + " -o " + outputDir;
  if (std::system(full_command.c_str()) != 0)
  {
    // in case of failure, try to run the python version
    auto hardcodedPythonPath = cbica::getExecutablePath() + "/OpenFederatedLearning/venv/bin/python"; // this needs to change for Windows (wonder what happens for macOS?)

    preparedataset_path = cbica::stringReplace(preparedataset_path, ".exe", "");
    preparedataset_path += ".py";
    full_command = hardcodedPythonPath + " " + preparedataset_path + " -i " + inputCSV + " -o " + outputDir;
    if (std::system(full_command.c_str()) != 0)
    {
      std::cerr << "There was an issue running PrepareDataset, contact 'software@cbica.upenn.edu' for troubleshooting.\n";
      return EXIT_FAILURE;
    }
  }

  // if it doesn't exit, run the fets_cli for inference
  auto fets_cli_path = getApplicationPath("FeTS_CLI");
  full_command = fets_cli_path + " -t 0 -d " + outputDir + "/DataForFeTS" + " -a " + archs + " -lF " + fusionMethod + " -g " + gpu_request;
  if (std::system(full_command.c_str()) != 0)
  {
    std::cerr << "There was an issue running FeTS_CLI, contact 'software@cbica.upenn.edu' for troubleshooting.\n";
    return EXIT_FAILURE;
  }

  std::cout << "Successfully finished.\n";

  return EXIT_SUCCESS;
}