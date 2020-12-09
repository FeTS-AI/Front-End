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

  auto hardcodedNativeModelWeightPath = getCaPTkDataDir() + "/fets";
  auto allArchs = cbica::subdirectoriesInDirectory(hardcodedNativeModelWeightPath);
  std::string allArchsString;
  for (size_t i = 0; i < allArchs.size(); i++)
  {
    allArchsString += allArchs[i] + ",";
  }
  allArchsString.pop_back();

  std::string dataDir, modelName, loggingDir, colName, archs = "3dresunet", fusionMethod = "STAPLE";

  parser.addRequiredParameter("d", "dataDir", cbica::Parameter::DIRECTORY, "Dir with Read/Write access", "Input data directory");
  parser.addRequiredParameter("t", "training", cbica::Parameter::BOOLEAN, "0 or 1", "Whether performing training or inference", "1==Train and 0==Inference");
  parser.addOptionalParameter("L", "LoggingDir", cbica::Parameter::DIRECTORY, "Dir with write access", "Location of logging directory");
  parser.addOptionalParameter("a", "archs", cbica::Parameter::STRING, allArchsString, "The architecture(s) to infer/train on", "Only a single architecture is supported for training", "Comma-separated values for multiple options", "Defaults to: " + archs);
  parser.addOptionalParameter("lF", "labelFuse", cbica::Parameter::STRING, "STAPLE,ITKVoting,SIMPLE,MajorityVoting", "The label fusion strategy to follow for multi-arch inference", "Comma-separated values for multiple options", "Defaults to: " + fusionMethod);
  parser.addOptionalParameter("g", "gpu", cbica::Parameter::BOOLEAN, "0-1", "Whether to run the process on GPU or not", "Defaults to '0'");
  parser.addOptionalParameter("c", "colName", cbica::Parameter::STRING, "", "Common name of collaborator", "Required for training");

  parser.addApplicationDescription("This is the CLI interface for FeTS");
  parser.addExampleUsage("-d /path/DataForFeTS -a deepMedic,nnUNet -lF STAPLE,ITKVoting,SIMPLE -g 1 -t 0", "This command performs inference using deepMedic,nnUNet using multiple fusion strategies on GPU and saves in data directory");
  parser.addExampleUsage("-d /path/DataForFeTS -t 1 -g 1 -c upenn", "This command starts training performs inference using deepMedic,nnUNet using multiple fusion strategies on GPU and saves in data directory");
  
  bool gpuRequested = false;
  bool trainingRequested = false;

  parser.getParameterValue("d", dataDir);
  parser.getParameterValue("t", trainingRequested);

  if (parser.isPresent("L"))
  {
    parser.getParameterValue("L", loggingDir);
  }
  else
  {
    loggingDir = dataDir + "/logs";
    cbica::createDir(loggingDir);
  }

  if (trainingRequested)
  {
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
  else
  {
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
  }
  if (parser.isPresent("g"))
  {
    parser.getParameterValue("g", gpuRequested);
  }

  std::string device_arg = " -md ";
  if (gpuRequested)
  {
    device_arg += "cuda";
  }
  else
  {
    device_arg += "cpu";
  }


  // convert everything to lower-case for easier comparison
  std::transform(archs.begin(), archs.end(), archs.begin(), ::tolower);
  std::transform(fusionMethod.begin(), fusionMethod.end(), fusionMethod.begin(), ::tolower);

  auto fetsApplicationPath = cbica::getExecutablePath();
  auto deepMedicExe = getApplicationPath("DeepMedic");

  auto archs_split = cbica::stringSplit(archs, ",");
  auto fusion_split = cbica::stringSplit(fusionMethod, ",");

  auto subjectDirs = cbica::subdirectoriesInDirectory(dataDir);

  if (trainingRequested && (archs_split.size() > 1))
  {
    std::cerr << "Training cannot be currently be performed on more than 1 architecture.\n";
    return EXIT_FAILURE;
  }

  std::string hardcodedPlanName,
    hardcodedOpenFLPath = fetsApplicationPath + "/OpenFederatedLearning/",
    hardcodedLabelFusionPath = fetsApplicationPath + "/LabelFusion/label_fusion",
    hardcodedModelWeightPath = hardcodedOpenFLPath + "/bin/federations/weights/", // start with the common location
    //hardcodedNativeModelWeightPath = hardcodedOpenFLPath + "/bin/federations/weights/native/", // the native weights are going in fets_data_dir/fets
    hardcodedPythonPath = hardcodedOpenFLPath + "/venv/bin/python"; // this needs to change for Windows (wonder what happens for macOS?)

  auto pythonEnvironmentFound = false;
  if (cbica::isFile(hardcodedPythonPath))
  {
    pythonEnvironmentFound = true;
  }

  if (!trainingRequested)
  {
    std::string subjectsWithMissingModalities, subjectsWithErrors; // string to store error cases
    
    std::cout << "Starting subject directory iteration...\n";
    for (size_t s = 0; s < subjectDirs.size(); s++) // iterate through all subjects
    {
      auto currentSubjectIsProblematic = false;
      std::string file_t1gd, file_t1, file_t2, file_flair;

      auto filesInDir = cbica::filesInDirectory(dataDir + "/" + subjectDirs[s]); // get files in current subject directory
      // iterate through all files and pick up individual modalities
      for (size_t f = 0; f < filesInDir.size(); f++)
      {
        if (file_t1gd.empty())
        {
          if ((filesInDir[f].find("_t1ce.nii.gz") != std::string::npos) || (filesInDir[f].find("_t1gd.nii.gz") != std::string::npos))
          {
            file_t1gd = filesInDir[f];
          }
        }
        if (file_t1.empty())
        {
          if (filesInDir[f].find("_t1.nii.gz") != std::string::npos)
          {
            file_t1 = filesInDir[f];
          }
        }
        if (file_t2.empty())
        {
          if (filesInDir[f].find("_t2.nii.gz") != std::string::npos)
          {
            file_t2 = filesInDir[f];
          }
        }
        if (file_flair.empty())
        {
          if ((filesInDir[f].find("_flair.nii.gz") != std::string::npos) || (filesInDir[f].find("_fl.nii.gz") != std::string::npos))
          {
            file_flair = filesInDir[f];
          }
        }
      }

      // ensure problematic cases are detected
      if (file_t1gd.empty())
      {
        subjectsWithMissingModalities += subjectDirs[s] + ",t1ce\n";
        currentSubjectIsProblematic = true;
      }
      if (file_t1.empty())
      {
        subjectsWithMissingModalities += subjectDirs[s] + ",t1\n";
        currentSubjectIsProblematic = true;
      }
      if (file_t2.empty())
      {
        subjectsWithMissingModalities += subjectDirs[s] + ",t2\n";
        currentSubjectIsProblematic = true;
      }
      if (file_flair.empty())
      {
        subjectsWithMissingModalities += subjectDirs[s] + ",flair\n";
        currentSubjectIsProblematic = true;
      }

      if (!currentSubjectIsProblematic) // proceed only if all modalities for the current subject are present
      {
        std::cout << "= Starting inference for subject: " << subjectDirs[s] << "\n";
        for (size_t a = 0; a < archs_split.size(); a++) // iterate through all requested architectures
        {
          if (archs_split[a] == "deepmedic") // special case 
          {
            std::cout << "== Starting inference using DeepMedic...\n";
            auto brainMaskFile = dataDir + "/" + subjectDirs[s] + "/" + subjectDirs[s] + "_deepmedic_seg.nii.gz";
            if (!cbica::isFile(brainMaskFile))
            {
              auto dm_tempOut = dataDir + "/" + subjectDirs[s] + "/dmOut/mask.nii.gz";
              auto fullCommand = deepMedicExe + " -md " + getCaPTkDataDir() + "/fets/deepMedic/saved_models/brainTumorSegmentation/ " +
                "-i " + file_t1 + "," +
                file_t1gd + "," +
                file_t2 + "," +
                file_flair + " -o " +
                dm_tempOut;

              if (std::system(fullCommand.c_str()) != 0)
              {
                std::cerr << "=== Couldn't complete the inference for deepmedic for subject " << subjectDirs[s] << ".\n";
                subjectsWithErrors += subjectDirs[s] + ",inference,deepmedic\n";
              }
              else
              {
                cbica::copyFile(dm_tempOut, brainMaskFile);
              }
            }
          } // deepmedic check
          else
          {
            auto fullCommandToRun = hardcodedPythonPath + " " + hardcodedOpenFLPath + "/bin/run_inference_from_flplan.py";
            auto args = " -d " + dataDir + device_arg + " -ld " + loggingDir + " -ip " + subjectDirs[s];
            if (pythonEnvironmentFound)
            {
              // check for all other models written in pytorch here

              // check between different architectures
              if (archs_split[a] == "3dunet")
              {
                // this is currently not defined
              }
              else if (archs_split[a] == "3dresunet")
              {
                std::cout << "== Starting inference using 3DResUNet...\n";
                hardcodedPlanName = "pt_3dresunet_brainmagebrats";
                auto hardcodedModelName = hardcodedPlanName + "_best.pbuf";
                auto allGood = true;
                if (!cbica::isFile((hardcodedModelWeightPath + "/" + hardcodedModelName))) // in case the "best" model is not present, use the "init" model that is distributed with FeTS installation
                {
                  auto hardcodedModelName = hardcodedPlanName + "_init.pbuf";
                  if (!cbica::isFile((hardcodedModelWeightPath + "/" + hardcodedModelName)))
                  {
                    std::cerr << "=== A compatible model weight file for the architecture '" << archs_split[a] << "' was not found. Please contact admin@fets.ai for help.\n";
                    allGood = false;
                  }
                }

                std::cout << "=== hardcodedModelName: " << hardcodedModelName << "\n";

                auto args_to_run = args + " -mwf " + hardcodedModelName
                  + " -p " + hardcodedPlanName + ".yaml";
                  //<< "-mwf" << hardcodedModelWeightPath // todo: doing customized solution above - change after model weights are using full paths for all
               
                if (std::system((fullCommandToRun + " " + args_to_run).c_str()) != 0)
                {
                  std::cerr << "=== Couldn't complete the inference for 3dresunet for subject " << subjectDirs[s] << ".\n";
                  subjectsWithErrors += subjectDirs[s] + ",inference,3dresunet\n";
                }
              } // end of 3dresunet check
              else
              {
                std::string hardcodedPlanName;
                if (archs_split[a].find("nnunet") != std::string::npos)
                {
                  hardcodedPlanName = "nnunet";
                  std::cout << "== Starting inference using nnUNet...\n";
                }
                else if (archs_split[a].find("deepscan") != std::string::npos)
                {
                  hardcodedPlanName = "deepscan";
                  std::cout << "== Starting inference using DeepScan...\n";
                }
                if (!hardcodedPlanName.empty())
                {
                  hardcodedPlanName += "_inference.yaml";
                  // structure according to what is needed - might need to create a function that can call run_inference_from_flplan for different hardcodedModelName
                  auto args_to_run = args + " -nmwf " + hardcodedNativeModelWeightPath // <abs path to folder containing all model weights folders> 
                    + " -p " + hardcodedPlanName
                    + " -pwai";

                  if (std::system((fullCommandToRun + " " + args_to_run).c_str()) != 0)
                  {
                    std::cerr << "=== Couldn't complete the inference for " << archs_split[a] << " for subject " << subjectDirs[s] << ".\n";
                    subjectsWithErrors += subjectDirs[s] + ",inference," + archs_split[a] + "\n";
                  }
                } // end of hardcodedPlanName check
              } // end of non-3dresunet check
            } // end of python check
          } // end of non-DM archs check
        } // end of archs_split

        /// fusion 
        if (pythonEnvironmentFound)
        {
          if (cbica::isFile(hardcodedLabelFusionPath))
          {
            std::cout << "== Starting label fusion...\n";
            auto filesInSubjectDir = cbica::filesInDirectory(dataDir + "/" + subjectDirs[s]);
            auto labelFusion_command = hardcodedPythonPath + " " + hardcodedLabelFusionPath + " ";
            std::string filesForFusion, dataForSegmentation = dataDir + "/" + subjectDirs[s] + "/SegmentationsForQC/";
            cbica::createDir(dataForSegmentation);

            for (size_t f = 0; f < filesInSubjectDir.size(); f++)
            {
              if (filesInSubjectDir[f].find("_seg.nii.gz") != std::string::npos) // find all files that have "_seg.nii.gz" in file name
              {
                if (filesInSubjectDir[f].find("final") == std::string::npos) // only do fusion for the files where "final" is not present
                {
                  auto fileToCopy = dataForSegmentation + cbica::getFilenameBase(filesInSubjectDir[f]) + ".nii.gz";
                  cbica::copyFile(filesInSubjectDir[f], fileToCopy);
                  filesForFusion += fileToCopy + ",";
                  std::remove(filesInSubjectDir[f].c_str());
                }
              }
            } // files loop in subject directory
            filesForFusion.pop_back(); // remove last ","

            for (size_t f = 0; f < fusion_split.size(); f++)
            {
              auto final_fused_file = dataForSegmentation + "/" + subjectDirs[s] + "_fused_" + fusion_split[f] + "_seg.nii.gz";
              auto full_fusion_command = labelFusion_command + "-inputs " + filesForFusion + " -classes 0,1,2,4 " // this needs to change after different segmentation algorithms are put in place
                + " -method " + fusion_split[f] + " -output " + final_fused_file;
              if (std::system(full_fusion_command.c_str()) != 0)
              {
                std::cerr << "=== Something went wrong with fusion for subject '" << subjectDirs[s] << "' using fusion method '" << fusion_split[f] << "'\n";
                subjectsWithErrors += subjectDirs[s] + ",fusion," + fusion_split[f] + "\n";
              }
            }
          } // end of label fusion script check
        } // end of python check
      } // end of currentSubjectIsProblematic 
    } // end of subjectDirs

    // provide error message
    if (!subjectsWithMissingModalities.empty())
    {
      std::cerr << "\nThe following subjects did not have all the 4 structural modalities to proceed with preprocessing:\nSubjectID,Modality\n" << subjectsWithMissingModalities;
    }
    if (!subjectsWithErrors.empty())
    {
      std::cerr << "\nThe following subjects were problematic:\nSubjectID,Application,Algorithm\n" << subjectsWithErrors;
    }
  } // end of trainingRequested check
  else // for training
  {
    std::string specialArgs, args, hardcodedModelName;
    if (trainingRequested)
    {
      specialArgs = "-col " + colName;
    }
    /// no longer checking for skull-stripping model - handled in PrepareDataset
    //if (modelName.find("_3dresunet_ss") != std::string::npos) // let's not worry about skull-stripping right now
    //{
    //  hardcodedPlanName = "pt_3dresunet_ss_brainmagebrats";
    //  hardcodedModelName = hardcodedModelWeightPath + hardcodedPlanName + "_best.pt"; // taken from https://github.com/FETS-AI/Models/blob/master/skullstripping/3dresunet/pt_3dresunet_ss_brainmagebrats_best.pt
    //  if (!trainingRequested)
    //  {
    //    specialArgs += "-nmwf " + hardcodedModelName;
    //  }
    //}
    //else
    {
      hardcodedPlanName = "pt_3dresunet_brainmagebrats"; // todo: this would need to changed based on the input arch in a future release
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

    args += device_arg;

    if (std::system((fullCommandToRun + " " + args + " " + specialArgs).c_str()) != 0)
    {
      std::cerr << "Couldn't complete the requested task.\n";
      return EXIT_FAILURE;
    }

  } // end of trainingRequested check

  std::cout << "Finished.\n";

  return EXIT_SUCCESS;
}


