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
        headers.push_back(cell);
      }
      else
      {
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


  return EXIT_SUCCESS;
}