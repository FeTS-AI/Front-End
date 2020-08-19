///////////////////////////////////////////////////////////////////////////////////////
// fTrainingPanel.h
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

#ifndef _fTrainingPanel_h_
#define _fTrainingPanel_h_


//#include "CAPTk.h"
#include "ui_fTrainingPanel.h"

#include "qsignalmapper.h"

/**
\class fTrainingPanel

\brief This class controls the elements in the Feature panel of the tab
*/
class fTrainingPanel : public QWidget, private Ui::fTrainingPanel
{
  Q_OBJECT

public:
  //! Constructor
  fTrainingPanel(QWidget * parent = 0);

  //! Destructor
  ~fTrainingPanel() {};
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

  //! Outputs which segmentation method has been selected and its corresponding config file (if any)
  std::pair< std::string, std::string > getSelectedSegmentationAlgorithms();

  std::string getOutputPath()
  {
    return m_outputPath.toStdString();
  }

  std::string getInputDirectoryPath()
  {
    return m_inputDataDirectory.toStdString();
  }

  std::string getCollaboratorName()
  {
    return m_collaboratorName->text().toStdString();
  }

  bool isGPUenabled()
  {
    return m_deviceControl->isChecked();
  }

signals:
  void m_btnComputeClicked();
  void helpClicked_SegmentationUsage(std::string);

public slots:
  void onComputeButtonClicked();
  void helpClicked();
  //void browseOutputFileName();
  void browseInputDirectory();
  void selectTrainingAlgorithm(int);
  void onCancelButtonClicked();
  void onOkButtonClicked();

private:
  QString mInputPathName, m_outputPath, m_inputDataDirectory;
  QString m_fetsDataDir = cbica::normalizePath(getCaPTkDataDir() + "/fets").c_str();
  std::map< int, std::map< int, std::string > > m_trainedModelConfigFile; //! variable to store the filename for associated model config to start re-training from
};


#endif
