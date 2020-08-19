#include "cbicaCmdParser.h"
#include "cbicaUtilities.h"

int main(int argc, char** argv)
{
  auto parser = cbica::CmdParser(argc, argv, "GreedyWrapper");
  parser.addRequiredParameter("f", "fixedImage", cbica::Parameter::FILE, "NIfTI", "The fixed image file");
  parser.addRequiredParameter("m", "movingImage", cbica::Parameter::FILE, "NIfTI", "The moving image file");
  parser.addRequiredParameter("o", "outputDir", cbica::Parameter::DIRECTORY, "None", "The output directory");
  parser.addOptionalParameter("d", "dimension", cbica::Parameter::INTEGER, "2-4", "The dimension of the images", "Default: 3");
  parser.addOptionalParameter("ov", "overWrite", cbica::Parameter::BOOLEAN, "0-1", "Overwrite previous results (if same parameters)", "Default: 0");

  std::string fixedImageFile, movingImageFile, outputDir;
  int dimension = 3;
  bool overWritePrevious = false;

  parser.getParameterValue("f", fixedImageFile);
  parser.getParameterValue("m", movingImageFile);
  parser.getParameterValue("o", outputDir);
  if (parser.isPresent("d"))
  {
    parser.getParameterValue("d", dimension);
  }

  if (!cbica::directoryExists(outputDir))
  {
    cbica::createDir(outputDir);
  }

  std::string common = "greedy" + 
#if WIN32
    std::string(".exe") +
#endif
    std::to_string(dimension);

  std::string commonFiles = "-n 100x50x10 -i " + fixedImageFile + " " + movingImageFile;
  std::string path, inputFile_base, ext;
  cbica::splitFileName(movingImageFile, path, inputFile_base, ext);
  
  /// change final parameter here
  std::string metric_command = "-m NCC 4x4x4";
  std::string currentMetric = "ncc444";

  std::string affine_out = outputDir + "/affine_noMask_" + currentMetric + ".mat";
  std::string deform_out = outputDir + "/deform_noMask_" + currentMetric + ".nii.gz";
  std::string deformInv_out = outputDir + "/deformInv_noMask_" + currentMetric + ".nii.gz";
  std::string finalOutput_affine = outputDir + "/" + inputFile_base + "_output_affine_pat2atl_noMask_" + currentMetric + ".nii.gz";
  std::string finalOutput = outputDir + "/" + inputFile_base + "_output_pat2atl_noMask_" + currentMetric + ".nii.gz";
  std::string finalOutput_Inv = outputDir + "/" + inputFile_base + "_output_atl2pat_noMask_" + currentMetric + ".nii.gz";

  auto command = "greedy.exe -d 3 -a " + metric_command + " " + common + " -ia-image-centers -o " + affine_out;
  if (!cbica::isFile(affine_out) || overWritePrevious)
  {
    std::cout << "[DEBUG] Started Affine ...\n";
    //std::cout << "[DEBUG] Command: " << command << "\n";
    std::system(command.c_str());
  }
  if (!cbica::isFile(deform_out) || overWritePrevious)
  {
    command = "greedy.exe -d 3 " + metric_command + " " + common + " -it " + affine_out + " -o " + deform_out + " -oinv " + deformInv_out;
    std::cout << "[DEBUG] Started Deformable ...\n";
    //std::cout << "[DEBUG] Command: " << command << "\n";
    std::system(command.c_str());
  }
  if (!cbica::isFile(finalOutput_affine) || overWritePrevious)
  {
    command = "greedy.exe -d 3 -rf " + fixedImageFile + " -rm " + movingImageFile+ " " + finalOutput_affine + " -ri LABEL 0.2vox -r " + affine_out;
    std::cout << "[DEBUG] Started Warping Affine ...\n";
    //std::cout << "[DEBUG] Command: " << command << "\n";
    std::system(command.c_str());
  }
  if (!cbica::isFile(finalOutput) || overWritePrevious)
  {
    command = "greedy.exe -d 3 -rf " + fixedImageFile + " -rm " + movingImageFile + " " + finalOutput + " -ri LABEL 0.2vox -r " + deform_out + " " + affine_out;
    std::cout << "[DEBUG] Started Warping ...\n";
    //std::cout << "[DEBUG] Command: " << command << "\n";
    std::system(command.c_str());
  }
  if (!cbica::isFile(finalOutput_Inv) || overWritePrevious)
  {
    command = "greedy.exe -d 3 -rf " + movingImageFile + " -rm " + fixedImageFile + " " + finalOutput_Inv + " -ri LABEL 0.2vox -r " + deformInv_out + " " + affine_out + ",-1";
    std::cout << "[DEBUG] Started WarpingInv ...\n";
    //std::cout << "[DEBUG] Command: " << command << "\n";
    std::system(command.c_str());
  }

  return EXIT_SUCCESS;
}