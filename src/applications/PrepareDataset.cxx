#include "cbicaCmdParser.h"
#include "cbicaUtilities.h"

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
        else if ((temp == "t2flair") || (temp == "flair") || (temp == "fl"))
        { 
          headers.push_back("FLAIR");
        }
      }
      else
      {
        if (headers.size() != 5)
        {
          std::cerr << "All required headers were not found in CSV. Please ensure the following are present: 'PatientID,T1,T1GD,T2,T2FLAIR'";
        }
        if (cell.find(" ") != std::string::npos)
        {
          std::cerr << "Please ensure that there are no spaces in the file paths.";
          return csvContents;
        }
        auto temp = cbica::stringReplace(cell, "\"", ""); // remove all double quotes
        temp = cbica::stringReplace(temp, "'", ""); // remove all single quotes
        currentRow[headers[j]] = temp;
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

  std::string bratsPipeline_exe = cbica::getExecutablePath() + "/BraTSPipeline";
#if WIN32
  bratsPipeline_exe += ".exe";
#endif
  if (!cbica::isFile(bratsPipeline_exe))
  {
    std::cerr << "BraTSPipeline was not found in the installation, cannot proceed.\n";
    return EXIT_FAILURE;
  }
  // iterate through all subjects
  for (size_t i = 0; i < csvContents.size(); i++)
  {
    auto subjectID = csvContents[i]["ID"];
    auto file_t1 = csvContents[i]["T1"];
    auto file_t1c = csvContents[i]["T1GD"];
    auto file_t2 = csvContents[i]["T2"];
    auto file_fl = csvContents[i]["FLAIR"];

    auto interimOutputDir = outputDir_qc + "/" + subjectID;
    auto finalSubjectOutputDir = outputDir_final + "/" + subjectID;
    cbica::createDir(finalSubjectOutputDir);

    auto command = bratsPipeline_exe + " -t1 " + csvContents[i]["T1"] + " -t1c " + csvContents[i]["T1GD"] + " -t2 " + csvContents[i]["T2"] + " -fl " + csvContents[i]["FLAIR"] + " -o " + interimOutputDir + " -s 1";
    if (std::system(command.c_str()) != 0)
    {
      std::cerr << "BraTSPipeline failed for subject " << csvContents[i]["ID"] << "\n";
    }
    else
    {
      cbica::copyFile(interimOutputDir + "/brain_T1CE.nii.gz", finalSubjectOutputDir + "/brain_t1gd.nii.gz");
      cbica::copyFile(interimOutputDir + "/brain_T1.nii.gz", finalSubjectOutputDir + "/brain_t1.nii.gz");
      cbica::copyFile(interimOutputDir + "/brain_T2.nii.gz", finalSubjectOutputDir + "/brain_t2.nii.gz");
      cbica::copyFile(interimOutputDir + "/brain_FLAIR.nii.gz", finalSubjectOutputDir + "/brain_flair.nii.gz");
    }
  }

  return EXIT_SUCCESS;
}