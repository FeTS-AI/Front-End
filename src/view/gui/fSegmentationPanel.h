///////////////////////////////////////////////////////////////////////////////////////
// fSegmentationPanel.h
//
// Copyright (c) 2018. All rights reserved.
// Section of Biomedical Image Analysis
// Center for Biomedical Image Computing and Analytics
// Department of Radiology
// Perelman School of Medicine
// University of Pennsylvania
//
// Contact details: software@cbica.upenn.edu
//
// License Agreement: https://www.med.upenn.edu/sbia/software-agreement.html
///////////////////////////////////////////////////////////////////////////////////////

#ifndef _fSegmentationPanel_h_
#define _fSegmentationPanel_h_


//#include "CAPTk.h"
#include "ui_fSegmentationPanel.h"

/**
\class fSegmentationPanel

\brief This class controls the elements in the Feature panel of the tab
*/
class fSegmentationPanel : public QWidget, private Ui::fSegmentationPanel
{
  Q_OBJECT

public:
  //! Constructor
  fSegmentationPanel(QWidget * parent = 0);

  //! Destructor
  ~fSegmentationPanel() {};
  // void  writeFeatureList(std::string Filename, std::vector< std::tuple<std::string, std::string, float>>featurevec);

  void setListner(void* lst)
  {
    //m_listener = lst;
  }

  void SetCurrentImagePath(const QString &inputPath)
  {
    mInputPathName = inputPath;
  }

  //! Sets the same temporary folder everywhere
  void setTempFolderLocation(const std::string& input_tempFolder)
  {
    m_tempFolderLocation = input_tempFolder;
  }

  //! Returns the data structure containing the segmentation check boxes
  std::vector< std::string > getSelectedSegmentationAlgorithms();

  //! Returns the data structure containing the fusion check boxes
  std::vector< std::string > getSelectedFusionAlgorithms();

  std::string getOutputPath()
  {
    return m_outputPath.toStdString();
  }

  std::string getInputDirectoryPath()
  {
    return m_inputDataDirectory.toStdString();
  }

  bool isGPUenabled()
  {
    return m_deviceControl->isChecked();
  }

  int getSelectedModel()
  {
    return m_availableModels->currentIndex();
  }

signals:
  void m_btnComputeClicked();
  void helpClicked_SegmentationUsage(std::string);

public slots:
  void onComputeButtonClicked();
  void helpClicked();
  //void browseOutputFileName();
  void browseInputDirectory();
  void selectSegmentationAlgorithm(int);
  void unSelectSegmentationAlgorithm(int);
  void onCancelButtonClicked();
  void onOkButtonClicked();
  //void callSegmentationAlgorithm

private:
  QString mInputPathName, m_outputPath, m_inputDataDirectory;
  QString m_fetsDataDir = this->GetFetsDataDir();

  int m_currentSegmentationAlgorithm; //! keeps a track of the selected segmentation algorithm

  std::map< int, std::map< int, std::string > > m_trainedModelConfigFile; //! variable to store the filename for associated trained model config
};


#endif
