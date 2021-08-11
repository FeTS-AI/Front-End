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
  auto hardcodedNativeModelWeightPath = fets_dataDir + "/fets";
  auto hardcodedFinalModelsWeightsPath = fets_dataDir + "/fets_consensus";

  std::string dataDir, outputDir, loggingDir, fusionMethod = "STAPLE", hardcodedPlanName = "fets_phase2_2";

  parser.addRequiredParameter("d", "dataDir", cbica::Parameter::DIRECTORY, "Dir with Read/Write access", "Input data directory");
  parser.addOptionalParameter("o", "outputDir", cbica::Parameter::DIRECTORY, "Dir with write access", "Location of logging directory");
  parser.addOptionalParameter("g", "gpu", cbica::Parameter::BOOLEAN, "0-1", "Whether to run the process on GPU or not", "Defaults to '0'");

  parser.addApplicationDescription("This is the CLI interface for FeTS");
  parser.addExampleUsage("-d /path/DataForFeTS -o /path/outputDir -g 1", "This command performs inference using the specific models and generates the output to send");
  
  bool gpuRequested = false;

  parser.getParameterValue("d", dataDir);
  parser.getParameterValue("o", outputDir);

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
    fullPlanPath = hardcodedOpenFLPath + "/bin/federations/plans/" + hardcodedPlanName + ".yaml"; // the full path to the plan that we want to use

  return EXIT_SUCCESS;
}