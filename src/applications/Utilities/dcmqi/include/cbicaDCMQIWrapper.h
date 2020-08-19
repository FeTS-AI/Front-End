#pragma once

namespace cbica
{
  void ConvertNiftiToDicomSeg(const std::string& niftiSegmentationFile, const std::string& dicomReferenceDirectory, const std::string& dicomSegMetaJSON, const std::string& outputDicomFile);
}