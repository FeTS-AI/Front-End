#include "cbicaCmdParser.h"
#include "cbicaLogging.h"
#include "cbicaITKSafeImageIO.h"
#include "cbicaUtilities.h"
#include "cbicaITKUtilities.h"
#include "CaPTkGUIUtils.h"

#include "BiasCorrection.hpp"

#include "itkMaskImageFilter.h"

#include <map>


int main(int argc, char** argv)
{
  cbica::CmdParser parser(argc, argv, "OpenFLCLI");

  
  std::cout << "Finished.\n";

  return EXIT_SUCCESS;
}


