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
    hardcodedOpenFLPlanPath = hardcodedOpenFLPath + "bin/federations/plans/fets_phase2_2.yaml",
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
  command_to_run = hardcodedPythonPath + " " + scriptToCall
    + " -WT " + hardcodedFinalModelsWeightsPath + "/overall -ET " + hardcodedFinalModelsWeightsPath + "/overall -TC " + hardcodedFinalModelsWeightsPath + "/overall "
    + "-pp " + hardcodedOpenFLPlanPath + " -op " + outputDir_overall + device_arg + " -dp " + dataDir + " -ptd";

  if (std::system(command_to_run.c_str()) != 0)
  {
    std::cerr << "The overall models did not run, please contact admin@fets.ai.\n\n";
    return EXIT_FAILURE;
  }

  auto outputDir_distinct = outputDir + "/distinct";
  cbica::createDir(outputDir_distinct);
  std::cout << "Starting distinct model scoring.\n";
  command_to_run = hardcodedPythonPath + " " + scriptToCall
    + " -WT " + hardcodedFinalModelsWeightsPath + "/WT -ET " + hardcodedFinalModelsWeightsPath + "/ET -TC " + hardcodedFinalModelsWeightsPath + "/TC "
    + "-pp " + hardcodedOpenFLPlanPath + " -op " + outputDir_overall + device_arg + " -dp " + dataDir + " -ptd";
  if (std::system(command_to_run.c_str()) == 0)
  {
    std::cerr << "The distinct models did not run, please contact admin@fets.ai.\n\n";
    return EXIT_FAILURE;
  }

  std::cout << "Starting inference for previous BraTS algorithms.\n";
  auto fetsCLIMainPath = getApplicationPath("FeTS_CLI");

  std::string fetsCLIArgs = " -d " + dataDir + " -t 0 -a deepMedic,nnUNet,DeepScan -g ";
  if (gpuRequested)
  {
    fetsCLIArgs += "1 ";
  }
  else
  {
    fetsCLIArgs += "0 ";
  }
  if (std::system((fetsCLIMainPath + fetsCLIArgs).c_str()) != 0)
  {
    std::cerr << "Previous BraTS algorithms did not run, hence validation results won't be generated for them.\n";
  }
  else
  {
    /// start validation of nnunet/deepscan/deepmedic on all validation cases
    auto validation_to_send = outputDir + "/validation_prev_brats.yaml",
      validation_internal = dataDir + "/validation_internal.yaml";

    if (!cbica::isFile(hardcodedPythonPath))
    {
      std::cerr << "The python virtual environment was not found, please refer to documentation to initialize it.\n";
      return EXIT_FAILURE;
    }

    auto regions_of_interest = { "WT", "TC", "ET" },
      measures_of_interest = { "Dice", "Hausdorff95", "Sensitivity", "Specificity" };

    auto yaml_config_to_send = YAML::Node();
    auto yaml_config_internal = YAML::Node();

    auto subjectsForConsideration = cbica::subdirectoriesInDirectory(dataDir);
    for (size_t i = 0; i < subjectsForConsideration.size(); i++)
    {
      auto currentSubjectDir = cbica::normPath(dataDir + "/" + subjectsForConsideration[i]);

      auto subject_id = subjectsForConsideration[i];
      auto subject_index_str = std::to_string(i);

      auto current_subject_folder = dataDir + "/" + subject_id;
      auto final_seg = current_subject_folder + "/" + subject_id + "_final_seg.nii.gz";
      std::map< std::string, std::string > archs_to_check;
      archs_to_check["deepmedic"] = current_subject_folder + "/SegmentationsForQC/" + subject_id + "_deepmedic_seg.nii.gz";
      archs_to_check["nnunet"] = current_subject_folder + "/SegmentationsForQC/" + subject_id + "_nnunet_seg.nii.gz";
      archs_to_check["deepscan"] = current_subject_folder + "/SegmentationsForQC/" + subject_id + "_deepscan_seg.nii.gz";
      if (!cbica::isFile(final_seg))
      {
        std::cerr << "The subject '" << subject_id << "' does not have a final_seg file present.\n";
      }
      else
      {
        using DefaultImageType = itk::Image< unsigned int, 3 >;
        auto final_seg_image = cbica::ReadImage< DefaultImageType >(final_seg);
        for (auto& current_arch : archs_to_check)
        {
          if (cbica::isFile(current_arch.second))
          {
            auto image_to_check = cbica::ReadImage< DefaultImageType >(current_arch.second);

            auto stats = cbica::GetBraTSLabelStatistics< DefaultImageType >(final_seg_image, image_to_check);

            for (auto& region : regions_of_interest)
            {
              for (auto& measure : measures_of_interest)
              {
                yaml_config_to_send[subject_index_str][current_arch.first][region][measure] = stats[region][measure];
                yaml_config_internal[subject_id][current_arch.first][region][measure] = stats[region][measure];
              } // end measure loop
            } // end region loop
          } // end file-check loop
        } // end arch-loop
      } // end final_seg check 
    } // end subject loop

    std::ofstream fout_int(validation_internal);
    fout_int << yaml_config_internal; // dump it back into the file
    fout_int.close();

    std::ofstream fout(validation_to_send);
    fout << yaml_config_to_send; // dump it back into the file
    fout.close();
  }

  std::cout << "FeTS Validation completed without errors. Please zip the following directory: '" << outputDir << "' and send a cloud storage link to admin@fets.ai.\n\n";

  return EXIT_SUCCESS;
}