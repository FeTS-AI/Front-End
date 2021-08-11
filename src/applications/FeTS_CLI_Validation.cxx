#include <map>

#include "cbicaCmdParser.h"
#include "cbicaLogging.h"
#include "cbicaITKSafeImageIO.h"
#include "cbicaUtilities.h"
#include "cbicaITKUtilities.h"
#include "CaPTkGUIUtils.h"

int main(int argc, char** argv)
{
  cbica::CmdParser parser(argc, argv, "FeTS_CLI");

  auto fets_dataDir = getCaPTkDataDir();
  auto hardcodedFinalModelsWeightsPath = fets_dataDir + "/fets_consensus";

  std::string dataDir, outputDir, loggingDir, fusionMethod = "STAPLE", hardcodedPlanName = "fets_phase2_2";

  parser.addRequiredParameter("d", "dataDir", cbica::Parameter::DIRECTORY, "Dir with Read/Write access", "Input data directory");
  parser.addOptionalParameter("o", "outputDir", cbica::Parameter::DIRECTORY, "Dir with write access", "Location of logging directory");
  parser.addOptionalParameter("g", "gpu", cbica::Parameter::BOOLEAN, "0-1", "Whether to run the process on GPU or not", "Defaults to '0'");

  parser.addApplicationDescription("This is the CLI interface for FeTS validation");
  parser.addExampleUsage("-d /path/DataForFeTS -o /path/outputDir -g 1", "This command performs inference using the specific models and generates the output to send");
  
  bool gpuRequested = false;

  parser.getParameterValue("d", dataDir);
  parser.getParameterValue("o", outputDir);
  cbica::createDir(outputDir);

  if (parser.isPresent("g"))
  {
    parser.getParameterValue("g", gpuRequested);
  }
  std::string device_arg = " -dev ";
  if (gpuRequested)
  {
    device_arg += "cuda";
  }
  else
  {
    device_arg += "cpu";
  }

  std::string
    fetsApplicationPath = cbica::getExecutablePath(),
    hardcodedOpenFLPath = fetsApplicationPath + "/OpenFederatedLearning/",
    hardcodedLabelFusionPath = fetsApplicationPath + "/LabelFusion/fusion_run",
    hardcodedModelWeightPath = hardcodedOpenFLPath + "/bin/federations/weights/", // start with the common location
    //hardcodedNativeModelWeightPath = hardcodedOpenFLPath + "/bin/federations/weights/native/", // the native weights are going in fets_data_dir/fets
    hardcodedPythonPath = hardcodedOpenFLPath + "/venv/bin/python", // this needs to change for Windows (wonder what happens for macOS?)
    hardcodedPythonPath_fusion = fetsApplicationPath + "/LabelFusion/venv/bin/python", // this needs to change for Windows (wonder what happens for macOS?)
    fullPlanPath = hardcodedOpenFLPath + "/bin/federations/plans/" + hardcodedPlanName + ".yaml", // the full path to the plan that we want to use
    scriptToCall = hardcodedOpenFLPath + "/submodules/fets_ai/Algorithms/fets/bin/brainmage_validation_scores_to_disk.py"; // the script that does the inference and scoring

  std::string command_to_run;
  
  auto outputDir_overall = outputDir + "/overall";
  cbica::createDir(outputDir_overall);
  std::cout << "Starting overall model scoring.\n";
  command_to_run = hardcodedPythonPath + scriptToCall
    + " -WT " + hardcodedFinalModelsWeightsPath + "/overall -ET " + hardcodedFinalModelsWeightsPath + "/overall -TC " + hardcodedFinalModelsWeightsPath + "/overall "
    + "-pp federations/plans/fets_phase2_2.yaml -op " + outputDir_overall + device_arg + " -dp " + dataDir + " -ptd";

  if (std::system(command_to_run.c_str()) != 0)
  {
    std::cerr << "The overall models did not run, please contact admin@fets.ai.\n\n";
    return EXIT_FAILURE;
  }

  auto outputDir_distinct = outputDir + "/distinct";
  cbica::createDir(outputDir_distinct);
  std::cout << "Starting distinct model scoring.\n";
  command_to_run = hardcodedPythonPath + scriptToCall
    + " -WT " + hardcodedFinalModelsWeightsPath + "/WT -ET " + hardcodedFinalModelsWeightsPath + "/ET -TC " + hardcodedFinalModelsWeightsPath + "/TC "
    + "-pp federations/plans/fets_phase2_2.yaml -op " + outputDir_overall + device_arg + " -dp " + dataDir + " -ptd";
  if (std::system(command_to_run.c_str()) == 0)
  {
    std::cerr << "The distinct models did not run, please contact admin@fets.ai.\n\n";
    return EXIT_FAILURE;
  }

  std::cout << "FeTS Validation completed without errors. Please zip the following directory: '" << outputDir << "' and send a cloud storage link to admin@fets.ai.\n\n";

  return EXIT_SUCCESS;
}