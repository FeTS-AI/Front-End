#include "cbicaCmdParser.h"
#include "cbicaLogging.h"
#include "cbicaITKSafeImageIO.h"
#include "cbicaUtilities.h"
#include "cbicaITKUtilities.h"
#include "CaPTkGUIUtils.h"

#include "itkMaskImageFilter.h"

#include <map>


int main(int argc, char** argv)
{
  cbica::CmdParser parser(argc, argv, "OpenFLCLI");

  auto allArchs = cbica::subdirectoriesInDirectory(getCaPTkDataDir() + "/fets");
  std::string allArchsString;
  for (size_t i = 0; i < allArchs.size(); i++)
  {
    allArchsString += allArchs[i] + ",";
  }
  allArchsString.pop_back();

  std::string dataDir, modelName, loggingDir, colName, archs, fusionOptions = "SIMPLE";

  parser.addRequiredParameter("d", "dataDir", cbica::Parameter::DIRECTORY, "Dir with Read/Write access", "Input data directory");
  parser.addRequiredParameter("m", "modelName", cbica::Parameter::FILE, "Model file", "Input model weights file");
  parser.addRequiredParameter("t", "training", cbica::Parameter::BOOLEAN, "0 or 1", "Whether performing training or inference", "1==Train and 0==Inference");
  parser.addRequiredParameter("L", "LoggingDir", cbica::Parameter::DIRECTORY, "Dir with write access", "Location of logging directory");
  parser.addRequiredParameter("a", "archs", cbica::Parameter::STRING, allArchsString, "The architecture(s) to infer/train on", "Only a single architecture is supported for training", "Comma-separated values for multiple options");
  parser.addOptionalParameter("lF", "labelFuse", cbica::Parameter::STRING, "SIMPLE,STAPLE,MajorityVoting", "The label fusion strategy to follow for multi-arch inference", "Comma-separated values for multiple options", "Defaults to: " + fusionOptions);
  parser.addOptionalParameter("g", "gpu", cbica::Parameter::BOOLEAN, "0-1", "Whether to run the process on GPU or not", "Defaults to '0'");
  parser.addOptionalParameter("c", "colName", cbica::Parameter::STRING, "", "Common name of collaborator", "Required for training");
  
  bool gpuRequested = false;
  bool trainingRequested = false;

  parser.getParameterValue("d", dataDir);
  parser.getParameterValue("m", modelName);
  parser.getParameterValue("L", loggingDir);
  parser.getParameterValue("a", archs);


  if (parser.isPresent("t"))
  {
    parser.getParameterValue("t", trainingRequested);
    if (parser.isPresent("c"))
    {
      parser.getParameterValue("c", colName);
    }
    else
    {
      std::cerr << "Collaborator name is required to beging training; please specify this using '-c'.\n";
      return EXIT_FAILURE;
    }
  }
  if (parser.isPresent("g"))
  {
    parser.getParameterValue("g", gpuRequested);
  }

  auto fetsApplicationPath = cbica::getExecutablePath();

  std::string hardcodedPlanName,
    hardcodedModelWeightPath = (fetsApplicationPath + "/OpenFederatedLearning/bin/federations/weights/"), // start with the common location
    hardcodedPythonPath = (fetsApplicationPath + "/OpenFederatedLearning/venv/bin/python"), // this needs to change for Windows (wonder what happens for macOS?)
    hardcodedModelName = "",
    args = "",
    specialArgs = "";

  if (trainingRequested)
  {
    specialArgs = "-col " + colName;
  }
  if (modelName.find("_3dresunet_ss") != std::string::npos)
  {
    hardcodedPlanName = "pt_3dresunet_ss_brainmagebrats";
    hardcodedModelName = hardcodedModelWeightPath + hardcodedPlanName + "_best.pt"; // taken from https://github.com/FETS-AI/Models/blob/master/skullstripping/3dresunet/pt_3dresunet_ss_brainmagebrats_best.pt
    if (!trainingRequested)
    {
      specialArgs += "-nmwf " + hardcodedModelName;
    }
  }
  else
  {
    hardcodedPlanName = "pt_3dresunet_brainmagebrats";
    auto hardcodedModelName = hardcodedPlanName + "_best.pbuf";
    if (!cbica::isFile((hardcodedModelWeightPath + "/" + hardcodedModelName)))
    {
      auto hardcodedModelName = hardcodedPlanName + "_init.pbuf";
      if (!cbica::isFile((hardcodedModelWeightPath + "/" + hardcodedModelName)))
      {
        std::cerr << "A compatible model weight file was not found. Please contact admin@fets.ai for help.\n";
        return EXIT_FAILURE;
      }
    }
    if (!trainingRequested)
    {
      specialArgs += "-mwf " + hardcodedModelName;
    }
  }

  // sanity checks
  //if (!cbica::isFile(hardcodedModelWeightPath.toStdString())) // todo: renable after model weights are using full paths for all
  //{
  //  ShowErrorMessage("The requested inference model was not found (it needs to be in ${FeTS_installDir}/bin/OpenFederatedLearning/bin/federations/weights/${planName}_best.pbuf");
  //  return;
  //}
  if (!cbica::isFile(hardcodedPythonPath))
  {
    std::cerr << "The python virtual environment was not found, please refer to documentation to initialize it.\n";
    return EXIT_FAILURE;
  }

  std::string fullCommandToRun = hardcodedPythonPath + " " + fetsApplicationPath;
  if (trainingRequested)
  {
    fullCommandToRun += "/OpenFederatedLearning/bin/run_inference_from_flplan.py";
  }
  else
  {
    fullCommandToRun += "/OpenFederatedLearning/bin/run_collaborator_from_flplan.py";
  }

  args += " -p " + hardcodedPlanName + ".yaml"
    //<< "-mwf" << hardcodedModelWeightPath // todo: doing customized solution above - change after model weights are using full paths for all
    + " -d " + dataDir
    + " -ld " + loggingDir;

  args += " -md ";
  if (gpuRequested)
  {
    args += "cuda";
  }
  else
  {
    args += "cpu";
  }

  if (std::system((fullCommandToRun + " " + args + " " + specialArgs).c_str()) != 0)
  {
    std::cerr << "Couldn't complete the requested task.\n";
    return EXIT_FAILURE;
  }
    
  std::cout << "Finished.\n";

  return EXIT_SUCCESS;
}


