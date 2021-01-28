//\*file  fMainWindow.cpp
//
//brief Implementation of fMainWindow class
//
//https://www.med.upenn.edu/sbia/software/ <br>
//software@cbica.upenn.edu
//
//Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
//See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html
//
//*/
/////
#include "fMainWindow.h"
#include "OutputInteractorStyleNavigator.h"
#include "SimpleImageManager.h"
#include "fHelpDialog.h"
//#include "EGFRvIIISurrogateIndex.h"
//#include "TrainingModule.h"
//#include "GeodesicSegmentation.h"
#include "N3BiasCorrection.h"
#include "SusanDenoising.h"
//#include "WhiteStripe.h"
//#include "PerfusionDerivatives.h"
//#include "PerfusionAlignment.h"
//#include "DiffusionDerivatives.h"
#include "ZScoreNormalizer.h"
//#include "PerfusionPCA.h"
//#include "SBRT_LungField.h"
//#include "SBRT_Nodule.h"
//#include "SBRT_Analysis.h"
#include "PreferencesDialog.h"

#include "cbicaITKSafeImageIO.h"
#include "itkFlipImageFilter.h"

#include "lddmm_common.h"
#include "lddmm_data.h"
#include "GreedyAPI.h"

#include <vnl/vnl_cost_function.h>
#include <vnl/vnl_random.h>
#include <vnl/algo/vnl_powell.h>
#include <vnl/algo/vnl_svd.h>
#include <vnl/vnl_trace.h>

#include "cbicaCmdParser.h"
#ifndef __APPLE__
//#include "LibraPreprocess.h"
#endif
#include "Registry.h"

//#include "DirectionalityEstimate.h"
#include "SlicerManager.h"
#include "Slicer.h"
#include "InteractorStyleNavigator.h"
#include "CaPTkUtils.h"
#include "CaPTkGUIUtils.h"

#include "vtkTransform.h"
#include "vtkImageMapToWindowLevelColors.h"
#include "vtkLookupTable.h"
#include "ComparisonViewerCommand.h"

#include "QtConcurrent/qtconcurrentrun.h"

#include "itkTranslationTransform.h"
#include "ApplicationPreferences.h"

#include "cbicaUtilities.h"

//#include "DicomSeriesReader.h"

// this function calls an external application from CaPTk in the most generic way while waiting for output
int fMainWindow::startExternalProcess(const QString &application, const QStringList &arguments)
{
  m_NumberOfUnfinishedExternalProcesses++;
  auto fullCommand = application.toStdString() + " " + arguments.join(" ").toStdString();
  cbica::Logging(loggerFile, fullCommand);
  int returnVal = std::system(fullCommand.c_str());
  m_NumberOfUnfinishedExternalProcesses--;
  return returnVal;

#ifdef _WIN32
  //QProcess process;
  //process.setStandardOutputFile((m_tempFolderLocation + "/process_" + application.toStdString() + ".log").c_str());
  //if (arguments.isEmpty())
  //{
  //  if (QFileInfo(application).completeSuffix() == "py")
  //  {
  //    process.start("python.exe " + application + ".py");
  //  }
  //  else
  //  {
  //    process.start(application);
  //  }
  //}
  //else
  //{
  //  if (QFileInfo(application).completeSuffix() == "py")
  //  {
  //    process.start("python.exe " + application + ".py", arguments);
  //  }
  //  else
  //  {
  //    process.start(application, arguments);
  //  }
  //}
  //process.write("exit\n\r");
  //process.waitForFinished(-1);
  //process.close();

  //return process.exitCode();
#else
  //std::string args_string = ""/*arguments.join(" ")*/, app_string = application.toStdString();
  //for (size_t i = 0; i < arguments.size(); i++)
  //{
  //  args_string += " " + arguments[i].toStdString();
  //}

  ////if (cbica::getFilenameExtension(application.toStdString()) == ".py")
  ////{
  ////  app_string = "python " + application.toStdString();
  ////}

  //return std::system((app_string + args_string).c_str());
#endif
}

int GetNumberOfDimensions(vtkImageData* input)
{
  int dim = 0;
#if VTK_MAJOR_VERSION <= 5
  int* extent = input->GetWholeExtent();
#else
  int* extent = input->GetExtent();
#endif
  if (extent[4] != extent[5])
  {
    dim = 3;
  }
  else if (extent[3] != extent[4])
  {
    dim = 2;
  }
  else if (extent[0] != extent[1])
  {
    dim = 1;
  }
  return dim;
}

inline std::string correctExtension(const std::string &inputFileName)
{
  std::string returnString = inputFileName, tempPath, tempBase, tempExt;
  cbica::splitFileName(returnString, tempPath, tempBase, tempExt);
  if (tempExt.empty())
  {
    returnString += NII_GZ_EXT;
  }

  return returnString;
}

fMainWindow::fMainWindow()
{

  setupUi(this);

  //! load preferences
  ApplicationPreferences::GetInstance()->DeSerializePreferences();
  ApplicationPreferences::GetInstance()->DisplayPreferences();

  //! comparison mode OFF at startup
  this->SetComparisonMode(false);

  this->bottomLayout = new QHBoxLayout();

  help_discussion = new QAction(this);
  help_forum = new QAction(this);
  help_bugs = new QAction(this);
  help_features = new QAction(this);
  help_download = new QAction(this);
  actionLoad_Recurrence_Images = new QAction(this);
  actionLoad_Nifti_Images = new QAction(this);
  actionLoad_Dicom_Images = new QAction(this);
  actionLoad_Nifti_ROI = new QAction(this);
  actionSave_Nifti_Images = new QAction(this);
  actionSave_Dicom_Images = new QAction(this);
  actionSave_ROI_Images = new QAction(this);
  actionSave_ROI_Dicom_Images = new QAction(this);
  actionExit = new QAction(this);
  actionAppEGFR = new QAction(this);
  actionAppRecurrence = new QAction(this);
  actionAppGeodesic = new QAction(this);
  actionAppGeodesicTraining = new QAction(this);
  actionHelp_Interactions = new QAction(this);
  actionAbout = new QAction(this);
  actionPreferences = new QAction(this);

  //---------------setting menu and status bar for the main window---------------
  this->setStatusBar(statusbar);

  menubar = new QMenuBar(this);
  menuFile = new QMenu("File");
  menuLoadFile = new QMenu("Load");
  menuSaveFile = new QMenu("Save");
  menuExit = new QMenu("Exit");
  menuLoadFileDicom = new QMenu("Dicom");
  menuLoadFileNifti = new QMenu("Nifti");
  menuFile->addMenu(menuLoadFile);
  menuFile->addMenu(menuSaveFile);
  menuApp = new QMenu("Applications");
  //menuDeepLearning = new QMenu("Deep Learning");
  menuPreprocessing = new QMenu("Preprocessing");
  menuHelp = new QMenu("Help");

  SaggitalViewWidget.reset(new QVTKOpenGLWidget(SaggitalWidget));
  AxialViewWidget.reset(new QVTKOpenGLWidget(AxialWidget));
  CoronalViewWidget.reset(new QVTKOpenGLWidget(CoronalWidget));

  SaggitalRenWin = vtkSmartPointer< vtkGenericOpenGLRenderWindow>::New();
  AxialRenWin = vtkSmartPointer< vtkGenericOpenGLRenderWindow>::New();
  CoronalRenWin = vtkSmartPointer< vtkGenericOpenGLRenderWindow>::New();

  SaggitalViewWidget->SetRenderWindow(SaggitalRenWin);
  AxialViewWidget->SetRenderWindow(AxialRenWin);
  CoronalViewWidget->SetRenderWindow(CoronalRenWin);

  SaggitalViewWidget->setMouseTracking(true);
  AxialViewWidget->setMouseTracking(true);
  CoronalViewWidget->setMouseTracking(true);

  SaggitalWidgetGridLayout->addWidget(SaggitalViewWidget.data(), 0, 0, 1, 1);
  AxialWidgetGridLayout->addWidget(AxialViewWidget.data(), 0, 0, 1, 1);
  CoronalWidgetGridLayout->addWidget(CoronalViewWidget.data(), 0, 0, 1, 1);

  QSizePolicy sizePolicy5(QSizePolicy::Preferred, QSizePolicy::Expanding);
  sizePolicy5.setHorizontalStretch(0);
  sizePolicy5.setVerticalStretch(0);

  preferenceDialog = new PreferencesDialog(nullptr);
  infoPanel = new fBottomImageInfoTip(centralwidget);
  imagesPanel = new fImagesPanel(); // New Images Panel
  m_tabWidget->addTab(imagesPanel, QString());
  //tumorPanel = new fTumorPanel();
  //m_tabWidget->addTab(tumorPanel, QString());
  drawingPanel = new fDrawingPanel();
  segmentationPanel = new fSegmentationPanel();
  trainingPanel = new fTrainingPanel();
  m_tabWidget->addTab(drawingPanel, QString());
  m_tabWidget->addTab(segmentationPanel, "Segmentation");
  m_tabWidget->addTab(trainingPanel, "Training");
  int minheight = drawingPanel->sizeHint().height() + 25;
  m_tabWidget->setMinimumHeight(minheight);
  m_tabWidget->setMaximumHeight(m_tabWidget->minimumHeight());

  m_toolTabdock->setWindowFlags(Qt::Window);

#ifdef Q_OS_WIN
  m_toolTabdock->setFeatures(QDockWidget::DockWidgetFloatable);
#else
  //TBD fix this - work around untill solved
  m_toolTabdock->setFeatures(QDockWidget::NoDockWidgetFeatures);
#endif
  m_toolTabdock->setWidget(m_tabWidget);
  overallGridLayout->addWidget(m_toolTabdock, 0, 0, 1, 3);

  QFrame * frame = new QFrame(this);
  sizePolicy5.setHeightForWidth(frame->sizePolicy().hasHeightForWidth());
  frame->setSizePolicy(sizePolicy5);
  frame->setFrameShape(QFrame::HLine);
  frame->setFrameShadow(QFrame::Sunken);

  overallGridLayout->addWidget(frame, 3, 0, 1, 3);

  this->setCentralWidget(centralwidget);
  AxialViewWidget->raise();
  CoronalViewWidget->raise();
  SaggitalViewWidget->raise();
  infoPanel->raise();
  m_tabWidget->raise();

  menuHelp->addAction(actionHelp_Interactions);
  //menuDownload = menuHelp->addMenu("Sample Data");
  //auto supportMenu = menuHelp->addMenu("Support Links");
  //menuHelp->addAction(actionAbout);

  //supportMenu->addAction(help_bugs);
  //supportMenu->addAction(help_download);

  menubar->addMenu(menuFile);
  menubar->addMenu(menuPreprocessing);
#ifndef PACKAGE_VIEWER
  //menubar->addMenu(menuApp);
#endif
  //menubar->addMenu(menuDeepLearning);
  menubar->addMenu(menuHelp);
  this->setMenuBar(menubar);

  menubar->addAction(menuFile->menuAction());
  menubar->addAction(menuPreprocessing->menuAction());
#ifndef PACKAGE_VIEWER
  //menubar->addAction(menuApp->menuAction());
#endif
  //menubar->addAction(menuDeepLearning->menuAction());
  menubar->addAction(menuHelp->menuAction());

  menuLoadFile->addAction(actionLoad_Nifti_Images);
  menuLoadFile->addAction(actionLoad_Nifti_ROI);
 // menuLoadFile->addAction(actionLoad_Dicom_Images);

  menuSaveFile->addAction(actionSave_Nifti_Images);
  menuSaveFile->addAction(actionSave_ROI_Images);

  menuFile->addAction(actionPreferences);
  menuFile->addAction(actionExit);

  //menuDownload->addAction("GreedyRegistration");
  m_tabWidget->setCurrentIndex(0);

  bottomLayout->addWidget(infoPanel);
  bottomLayout->addStretch();
  bottomLayout->addWidget(preferencesGroupBox);
  overallGridLayout->addLayout(bottomLayout, 4, 0, 2, 3);

  std::string nonNativeAppPaths_wrap = std::string(CAPTK_APP_LIST_PY_GUI);
  if (nonNativeAppPaths_wrap[0] == ' ')
  {
    nonNativeAppPaths_wrap.erase(0, 1);
  }
  //nonNativeAppPaths_wrap = nonNativeAppPaths_wrap + " itksnap";
  m_pyGUIApps = cbica::stringSplit(nonNativeAppPaths_wrap, " ");
  nonNativeAppPaths_wrap = std::string(CAPTK_APP_LIST_PY_CLI);
  if (nonNativeAppPaths_wrap[0] == ' ')
  {
    nonNativeAppPaths_wrap.erase(0, 1);
  }

  m_pyCLIApps = cbica::stringSplit(nonNativeAppPaths_wrap, " ");
  size_t allAppCounter = 0;
  for (size_t i = 0; i < m_pyGUIApps.size(); i++)
  {
    if (m_pyGUIApps[i] == "confetti")
    {
      m_pyGUIApps[i] = "ConfettiGUI";
    }
    if ((m_pyGUIApps[i] == "librabatch") || (m_pyGUIApps[i] == "librasingle"))
    {
      m_pyGUIApps[i] = "libra";
    }
    if (m_pyGUIApps[i] == "SBRT_Segment")
    {
      m_pyGUIApps[i] = "SBRT_Lung_Segment";
    }
    if (m_pyGUIApps[i] == "SBRT_Analyze")
    {
      m_pyGUIApps[i] = "SBRT_Lung_Analyze";
    }

    m_allNonNativeApps[m_pyGUIApps[i]] = getApplicationPath(m_pyGUIApps[i]);
    allAppCounter++;
  }

  for (size_t i = 0; i < m_pyCLIApps.size(); i++)
  {
    m_allNonNativeApps[m_pyCLIApps[i]] = getApplicationPath(m_pyCLIApps[i]);
  }

  // TBD: this needs to be controlled from CMake and not hard-coded here
  std::string brainAppList = "";// " EGFRvIIISVMIndex EGFRvIIISurrogateIndex RecurrenceEstimator PseudoProgressionEstimator SurvivalPredictor MolecularSubtypePredictor PopulationAtlases WhiteStripe confetti";
  std::string breastAppList = "";

#ifndef __APPLE__
  //breastAppList = " librasingle librabatch breastSegment texturePipeline";
#endif

  std::string lungAppList = "";// " LungField Nodule Analysis";
  //std::string miscAppList = " DirectionalityEstimate DiffusionDerivatives PerfusionAlignment PerfusionDerivatives PerfusionPCA TrainingModule";
  std::string miscAppList = "";// " DirectionalityEstimate DiffusionDerivatives TrainingModule";
  std::string segAppList = "";// " itksnap GeodesicSegmentation GeodesicTrainingSegmentation deepmedic_tumor deepmedic_brain";
  std::string preProcessingAlgos = " BiasCorrect-N3 Denoise-SUSAN GreedyRegistration HistogramMatching ZScoringNormalizer deepmedic_brain BraTSPipeline";
  //#ifndef __APPLE__
//  preProcessingAlgos += " breastNormalize";
//#endif
  //std::string deepLearningAlgos = " deepmedic_tumor deepmedic_brain";

  if (!brainAppList.empty())
  {
    vectorOfGBMApps = populateStringListInMenu(brainAppList, this, menuApp, "Glioblastoma", false);
    menuApp->addSeparator();
  }
  if (!breastAppList.empty())
  {
    vectorOfBreastApps = populateStringListInMenu(breastAppList, this, menuApp, "Breast Cancer", false);
    menuApp->addSeparator();
  }
  if (!lungAppList.empty())
  {
    vectorOfLungApps = populateStringListInMenu(lungAppList, this, menuApp, "Lung Cancer", false);
    menuApp->addSeparator();
  }
  if (!segAppList.empty())
  {
    vectorOfSegmentationApps = populateStringListInMenu(segAppList, this, menuApp, "Segmentation", false);
  }
  if (!miscAppList.empty())
  {
    vectorOfMiscApps = populateStringListInMenu(miscAppList, this, menuApp, "Miscellaneous", false);
  }
  if (!preProcessingAlgos.empty())
  {
    vectorOfPreprocessingActionsAndNames = populateStringListInMenu(preProcessingAlgos, this, menuPreprocessing, "", false);
  }
  //if (!deepLearningAlgos.empty())
  //{
  //  vectorOfDeepLearningActionsAndNames = populateStringListInMenu(deepLearningAlgos, this, menuDeepLearning, "", false);

  //  auto temp = populateStringListInMenu(" ", this, menuDeepLearning, "Breast", false);
  //  temp = populateStringListInMenu(" ", this, menuDeepLearning, "Lung", false);
  //  menuDeepLearning->addSeparator();
  //  temp = populateStringListInMenu(" ", this, menuDeepLearning, "Training", false);
  //}

  //menuDownload->addAction("All");
  //for (const auto &currentActionAndName : vectorOfGBMApps)
  //{
  //  if (currentActionAndName.name != "Glioblastoma")
  //  {
  //    if (currentActionAndName.name == "confetti")
  //    {
  //      menuDownload->addAction("Confetti");
  //    }
  //    else
  //    {
  //      menuDownload->addAction(currentActionAndName.name.c_str());
  //    }
  //  }
  //}

  bool libraCheck = false;
  for (const auto &currentActionAndName : vectorOfBreastApps)
  {
    if (!libraCheck)
    {
      if (currentActionAndName.name.find("libra") != std::string::npos)
      {
        libraCheck = true;
        //menuDownload->addAction("LIBRA");
      }
    }
    //if (currentActionAndName.name != "Breast Cancer")
    //{
    //  if (!libraCheck)
    //  {
    //    if (currentActionAndName.name.find("libra") != std::string::npos)
    //    {
    //      libraCheck = true;
    //      menuDownload->addAction("LIBRA");
    //    }
    //  }
    //}
  }

  bool sbrtCheck = false;
  for (const auto &currentActionAndName : vectorOfLungApps)
  {
    if (currentActionAndName.name != "Lung Cancer")
    {
      if (!sbrtCheck)
      {
        if ((currentActionAndName.name.find("LungField") != std::string::npos) ||
          (currentActionAndName.name.find("Nodule") != std::string::npos) ||
          (currentActionAndName.name.find("Analysis") != std::string::npos))
        {
          sbrtCheck = true;
          //menuDownload->addAction("LungCancer");
        }
      }
    }
  }

  for (const auto &currentActionAndName : vectorOfMiscApps)
  {
    if (currentActionAndName.name != "Miscellaneous")
    {
      if ((currentActionAndName.name != "itksnap") && (currentActionAndName.name != "deepmedic"))
      {
        if (!currentActionAndName.name.empty())
        {
          //menuDownload->addAction(currentActionAndName.name.c_str());
        }
      }
    }
  }


  m_imagesTable = imagesPanel->GetImagesTable();
  m_nonVisImagesTable = imagesPanel->GetNonViewingImagesTable();
  assert(m_imagesTable != NULL);
  assert(m_nonVisImagesTable != NULL);

  mSequenceNumber = 0;

  t1cePath = "";

  QString msg = tr(EXE_NAME);
#ifdef SW_VER
  msg += " - v" + tr(SW_VER);
#endif
  this->setWindowTitle(msg);

#ifdef Q_OS_WIN32
  currentPlatform = "windows";
#endif

  mInputPathName = "";
  //mMainWidget = this;
  mCurrentSelectedImageId = "";
  mCurrentPickedImageId = "";
  mCurrentPickedImageIndex = 0;
  mCurrentNearPoints = 0;
  mCurrentFarPoints = 0;
  mCurrentInitPoints = 0;
  mCustomImageToThreshold = itk::Image< short, 3 >::New();
  mProjectVariant = std::string(PROJECT_VARIANT);

  connect(&registrationPanel,
    SIGNAL(RegistrationSignal(std::string, std::vector<std::string>, std::vector<std::string>, std::vector<std::string>, std::string, bool, bool, bool, std::string, std::string)),
    this,
    SLOT(Registration(std::string, std::vector<std::string>, std::vector<std::string>, std::vector<std::string>, std::string, bool, bool, bool, std::string, std::string)));

  cbica::createDir(loggerFolder);
  m_tempFolderLocation = loggerFolder + "tmp_" + cbica::getCurrentProcessID();
  if (cbica::directoryExists(m_tempFolderLocation))
  {
    auto temp = cbica::stringSplit(cbica::getCurrentLocalTime(), ":");
    m_tempFolderLocation += temp[0] + temp[1] + temp[2] + "/";
  }
  cbica::createDir(m_tempFolderLocation);

  mLandmarks = new Landmarks(LANDMARK_TYPE::DEFAULT);
  mSeedPoints = new Landmarks(LANDMARK_TYPE::TUMOR_POINTS);
  mTissuePoints = new Landmarks(LANDMARK_TYPE::TISSUE_POINTS);

  mMask = vtkSmartPointer<vtkImageData>::New();

  mSlicerManagers.resize(0);

  image4DSlider->setEnabled(false);
  image4DSlider->setValue(0);

  connect(imagesPanel, SIGNAL(sigOverlayCheckBoxChanged(int)), this, SLOT(overlayUseStateChanged(int)));
  connect(imagesPanel, SIGNAL(sigOverlaySliderChanged(int)), this, SLOT(overlaySliderChanged(int)));
  connect(imagesPanel, SIGNAL(sigOverlayChanged()), this, SLOT(overlayChanged()));
  connect(imagesPanel, SIGNAL(sigTheiaClicked()), this, SLOT(ApplicationTheia()));
  connect(imagesPanel, SIGNAL(CompareModeToggled(bool)), this, SLOT(EnableComparisonMode(bool)));
  connect(imagesPanel, SIGNAL(sigImageModalityChanged(int)), this, SLOT(imageModalityChanged(int)));
  connect(imagesPanel, SIGNAL(helpClicked_Interaction(std::string)), this, SLOT(help_contextual(std::string)));

  connect(image4DSlider, SIGNAL(valueChanged(int)), this, SLOT(imageSliderChanged()));
  SetPresetComboBox();

  // init the sliders
  verticalSliders.push_back(AxialViewSlider);
  verticalSliders.push_back(CoronalViewSlider);
  verticalSliders.push_back(SaggitalViewSlider);
  for (int i = 0; i < 3; i++)
  {
    verticalSliders[i]->hide();
  }
  connect(AxialViewSlider, SIGNAL(valueChanged(int)), this, SLOT(AxialViewSliderChanged()));
  connect(CoronalViewSlider, SIGNAL(valueChanged(int)), this, SLOT(CoronalViewSliderChanged()));
  connect(SaggitalViewSlider, SIGNAL(valueChanged(int)), this, SLOT(SaggitalViewSliderChanged()));

  connect(actionLoad_Recurrence_Images, SIGNAL(triggered()), this, SLOT(openImages()));
  connect(actionLoad_Nifti_Images, SIGNAL(triggered()), this, SLOT(openImages()));
  //connect(actionLoad_Dicom_Images, SIGNAL(triggered()), this, SLOT(openDicomImages()));

  connect(actionSave_ROI_Images, SIGNAL(triggered()), this, SLOT(SaveDrawing()));
  connect(actionSave_ROI_Dicom_Images, SIGNAL(triggered()), this, SLOT(SaveDicomDrawing()));
  connect(actionSave_Nifti_Images, SIGNAL(triggered()), this, SLOT(SaveImage()));
  connect(actionSave_Dicom_Images, SIGNAL(triggered()), this, SLOT(SaveDicomImage()));
  connect(actionPreferences, SIGNAL(triggered()), this, SLOT(OnPreferencesMenuClicked()));

  connect(actionLoad_Nifti_ROI, SIGNAL(triggered()), this, SLOT(LoadDrawing()));

  connect(actionExit, SIGNAL(triggered()), this, SLOT(close()));
  connect(actionAbout, SIGNAL(triggered()), this, SLOT(about()));
  connect(actionHelp_Interactions, SIGNAL(triggered()), this, SLOT(help_Interactions()));
  connect(help_bugs, SIGNAL(triggered()), this, SLOT(help_BugTracker()));

  //connect(menuDownload, SIGNAL(triggered(QAction*)), this, SLOT(help_Download(QAction*)));

  connect(&mHelpTutorial, SIGNAL(skipTutorialOnNextRun(bool)), this, SLOT(skipTutorial(bool)));

  for (size_t i = 0; i < vectorOfGBMApps.size(); i++)
  {
    if (vectorOfGBMApps[i].name.find("EGFRvIIISurrogate") != std::string::npos)
    {
      vectorOfGBMApps[i].action->setText("  Glioblastoma EGFRvIII Surrogate Index (PHI Estimator)"); // TBD set at source
      connect(vectorOfGBMApps[i].action, SIGNAL(triggered()), this, SLOT(ApplicationEGFR()));
    }
    if (vectorOfGBMApps[i].name.find("EGFRvIIISVM") != std::string::npos)
    {
      vectorOfGBMApps[i].action->setText("  Glioblastoma EGFRvIII SVM Index"); // TBD set at source
      connect(vectorOfGBMApps[i].action, SIGNAL(triggered()), this, SLOT(ApplicationEGFRvIIISVM()));
    }
    else if (vectorOfGBMApps[i].name.find("Recurrence") != std::string::npos)
    {
      vectorOfGBMApps[i].action->setText("  Glioblastoma Infiltration Index"); // TBD set at source
      connect(vectorOfGBMApps[i].action, SIGNAL(triggered()), this, SLOT(ApplicationRecurrence()));
    }
    else if (vectorOfGBMApps[i].name.find("PseudoProgression") != std::string::npos)
    {
      vectorOfGBMApps[i].action->setText("  Glioblastoma PseudoProgression Infiltration Index"); // TBD set at source
      connect(vectorOfGBMApps[i].action, SIGNAL(triggered()), this, SLOT(ApplicationPseudoProgression()));
    }
    else if (vectorOfGBMApps[i].name.find("Survival") != std::string::npos)
    {
      vectorOfGBMApps[i].action->setText("  Glioblastoma Survival Prediction Index"); // TBD set at source
      connect(vectorOfGBMApps[i].action, SIGNAL(triggered()), this, SLOT(ApplicationSurvival()));
    }
    else if (vectorOfGBMApps[i].name.find("PopulationAtlases") != std::string::npos)
    {
      vectorOfGBMApps[i].action->setText("  Population Atlas"); // TBD set at source
      connect(vectorOfGBMApps[i].action, SIGNAL(triggered()), this, SLOT(ApplicationPopulationAtlas()));
    }
    else if (vectorOfGBMApps[i].name.find("ImagingSubtype") != std::string::npos)
    {
      vectorOfGBMApps[i].action->setText("  Glioblastoma Imaging Subtype Predictor"); // TBD set at source
      connect(vectorOfGBMApps[i].action, SIGNAL(triggered()), this, SLOT(ApplicationImagingSubtype()));
    }
    else if (vectorOfGBMApps[i].name.find("MolecularSubtypePredictor") != std::string::npos)
    {
      vectorOfGBMApps[i].action->setText("  Glioblastoma Molecular Subtype Predictor"); // TBD set at source
      connect(vectorOfGBMApps[i].action, SIGNAL(triggered()), this, SLOT(ApplicationMolecularSubtype()));
    }
    else if (vectorOfGBMApps[i].name.find("WhiteStripe") != std::string::npos)
    {
      vectorOfGBMApps[i].action->setText("  WhiteStripe Normalization"); // TBD set at source
      connect(vectorOfGBMApps[i].action, SIGNAL(triggered()), this, SLOT(ApplicationWhiteStripe()));
    }
    else if (vectorOfGBMApps[i].name.find("confetti") != std::string::npos)
    {
      vectorOfGBMApps[i].action->setText("  Confetti"); //TBD set at source
      connect(vectorOfGBMApps[i].action, SIGNAL(triggered()), this, SLOT(ApplicationConfetti()));
    }
    else if (vectorOfGBMApps[i].name.find("DirectionalityEstimate") != std::string::npos)
    {
      vectorOfGBMApps[i].action->setText("  Directionality Estimator"); //TBD set at source
      connect(vectorOfGBMApps[i].action, SIGNAL(triggered()), this, SLOT(ApplicationDirectionality()));
    }
  }

  for (size_t i = 0; i < vectorOfBreastApps.size(); i++)
  {
    if (vectorOfBreastApps[i].name.find("librasingle") != std::string::npos)
    {
      vectorOfBreastApps[i].action->setText("  Breast Density Estimator (LIBRA) SingleImage"); //TBD set at source
      connect(vectorOfBreastApps[i].action, SIGNAL(triggered()), this, SLOT(ApplicationLIBRASingle()));
    }
    else if (vectorOfBreastApps[i].name.find("librabatch") != std::string::npos)
    {
      vectorOfBreastApps[i].action->setText("  Breast Density Estimator (LIBRA) BatchMode"); //TBD set at source
      connect(vectorOfBreastApps[i].action, SIGNAL(triggered()), this, SLOT(ApplicationLIBRABatch()));
    }
    else if (vectorOfBreastApps[i].name.find("breastSegment") != std::string::npos)
    {
      vectorOfBreastApps[i].action->setText("  Breast Segmentation"); //TBD set at source
      connect(vectorOfBreastApps[i].action, SIGNAL(triggered()), this, SLOT(ApplicationBreastSegmentation()));
    }
    else if (vectorOfBreastApps[i].name.find("texturePipeline") != std::string::npos)
    {
      vectorOfBreastApps[i].action->setText("  Texture Feature Pipeline"); //TBD set at source
      connect(vectorOfBreastApps[i].action, SIGNAL(triggered()), this, SLOT(ApplicationTexturePipeline()));
    } 
  }

  for (size_t i = 0; i < vectorOfLungApps.size(); i++)
  {
    if (vectorOfLungApps[i].name.find("LungField") != std::string::npos)
    {
      vectorOfLungApps[i].action->setText("  Lung Field Segmentation");
      connect(vectorOfLungApps[i].action, SIGNAL(triggered()), this, SLOT(ApplicationSBRTLungField()));
    }
    if (vectorOfLungApps[i].name.find("Nodule") != std::string::npos)
    {
      vectorOfLungApps[i].action->setText("  Lung Nodule Segmentation");
      connect(vectorOfLungApps[i].action, SIGNAL(triggered()), this, SLOT(ApplicationSBRTNodule()));
    }
    if (vectorOfLungApps[i].name.find("Analysis") != std::string::npos)
    {
      vectorOfLungApps[i].action->setText("  Prognostic Modeling");
      connect(vectorOfLungApps[i].action, SIGNAL(triggered()), this, SLOT(ApplicationSBRTAnalysis()));
    }
  }

  for (size_t i = 0; i < vectorOfSegmentationApps.size(); i++)
  {
    if (vectorOfSegmentationApps[i].name.find("itksnap") != std::string::npos)
    {
      vectorOfSegmentationApps[i].action->setText("  ITK-SNAP"); //TBD set at source
      connect(vectorOfSegmentationApps[i].action, SIGNAL(triggered()), this, SLOT(ApplicationITKSNAP()));
    }
    else if (vectorOfSegmentationApps[i].name.find("GeodesicSegmentation") != std::string::npos)
    {
      vectorOfSegmentationApps[i].action->setText("  Geodesic Segmentation"); // TBD set at source
      connect(vectorOfSegmentationApps[i].action, SIGNAL(triggered()), this, SLOT(ApplicationGeodesic()));
    }
    else if (vectorOfSegmentationApps[i].name.find("GeodesicTrainingSegmentation") != std::string::npos)
    {
      vectorOfSegmentationApps[i].action->setText("  Geodesic Training Segmentation"); // TBD set at source
      connect(vectorOfSegmentationApps[i].action, SIGNAL(triggered()), this, SLOT(ApplicationGeodesicTraining()));
    }
    else if (vectorOfSegmentationApps[i].name.find("deepmedic_tumor") != std::string::npos)
    {
      vectorOfSegmentationApps[i].action->setText("  Brain Tumor Segmentation (DeepLearning)"); // TBD set at source
      connect(vectorOfSegmentationApps[i].action, &QAction::triggered, this, [this] { ApplicationDeepMedicSegmentation(fDeepMedicDialog::Tumor); });
    }
    else if (vectorOfSegmentationApps[i].name.find("deepmedic_brain") != std::string::npos)
    {
      vectorOfSegmentationApps[i].action->setText("  Skull Stripping (DeepLearning)"); // TBD set at source
      connect(vectorOfSegmentationApps[i].action, &QAction::triggered, this, [this] { ApplicationDeepMedicSegmentation(fDeepMedicDialog::SkullStripping); });
    }
  }

  for (size_t i = 0; i < vectorOfMiscApps.size(); i++)
  {
    if (vectorOfMiscApps[i].name.find("DirectionalityEstimate") != std::string::npos)
    {
      vectorOfMiscApps[i].action->setText("  Directionality Estimator"); //TBD set at source
      connect(vectorOfMiscApps[i].action, SIGNAL(triggered()), this, SLOT(ApplicationDirectionality()));
    }
    else if (vectorOfMiscApps[i].name.find("PerfusionPCA") != std::string::npos)
    {
      vectorOfMiscApps[i].action->setText("  Perfusion PCA"); //TBD set at source
      connect(vectorOfMiscApps[i].action, SIGNAL(triggered()), this, SLOT(ApplicationPCA()));
    }
    else if (vectorOfMiscApps[i].name.find("PerfusionDerivatives") != std::string::npos)
    {
      vectorOfMiscApps[i].action->setText("  Perfusion Derivatives"); //TBD set at source
      connect(vectorOfMiscApps[i].action, SIGNAL(triggered()), this, SLOT(ApplicationPerfusionMeasuresCalculation()));
    }
    else if (vectorOfMiscApps[i].name.find("PerfusionAlignment") != std::string::npos)
    {
      vectorOfMiscApps[i].action->setText("  Perfusion Alignment"); //TBD set at source
      connect(vectorOfMiscApps[i].action, SIGNAL(triggered()), this, SLOT(ApplicationPerfusionAlignmentCalculation()));
    }
    else if (vectorOfMiscApps[i].name.find("DiffusionDerivatives") != std::string::npos)
    {
      vectorOfMiscApps[i].action->setText("  Diffusion Derivatives"); //TBD set at source
      connect(vectorOfMiscApps[i].action, SIGNAL(triggered()), this, SLOT(ApplicationDiffusionMeasuresCalculation()));
    }
    else if (vectorOfMiscApps[i].name.find("TrainingModule") != std::string::npos)
    {
      vectorOfMiscApps[i].action->setText("  Training Module"); //TBD set at source
      connect(vectorOfMiscApps[i].action, SIGNAL(triggered()), this, SLOT(ApplicationTrainingModule()));
    }
  }

  // add a single function for all preprocessing steps, this function will check for the specific names and then initiate that algorithm
  for (size_t i = 0; i < vectorOfPreprocessingActionsAndNames.size(); i++)
  {
    if (vectorOfPreprocessingActionsAndNames[i].name.find("Denoise") != std::string::npos)
    {
      connect(vectorOfPreprocessingActionsAndNames[i].action, SIGNAL(triggered()), this, SLOT(ImageDenoising()));
    }

    else if (vectorOfPreprocessingActionsAndNames[i].name.find("BiasCorrect") != std::string::npos)
    {
      vectorOfPreprocessingActionsAndNames[i].action->setText("BiasCorrection");//TBD set at source
      connect(vectorOfPreprocessingActionsAndNames[i].action, SIGNAL(triggered()), this, SLOT(ImageBiasCorrection()));
    }
    else if (vectorOfPreprocessingActionsAndNames[i].name.find("GreedyRegistration") != std::string::npos)
    {
      vectorOfPreprocessingActionsAndNames[i].action->setText("Registration");
      connect(vectorOfPreprocessingActionsAndNames[i].action, SIGNAL(triggered()), this, SLOT(ImageRegistration()));
    }
    else if (vectorOfPreprocessingActionsAndNames[i].name.find("HistogramMatching") != std::string::npos)
    {
      connect(vectorOfPreprocessingActionsAndNames[i].action, SIGNAL(triggered()), this, SLOT(ImageHistogramMatching()));
    }
    else if (vectorOfPreprocessingActionsAndNames[i].name.find("DeepMedicNormalizer") != std::string::npos)
    {
      vectorOfPreprocessingActionsAndNames[i].action->setText("Z-Scoring Normalizer"); // TBD set at source
      connect(vectorOfPreprocessingActionsAndNames[i].action, SIGNAL(triggered()), this, SLOT(ImageDeepMedicNormalizer()));
    }
    else if (vectorOfPreprocessingActionsAndNames[i].name.find("SkullStripping") != std::string::npos)
    {
      vectorOfPreprocessingActionsAndNames[i].action->setText("Skull Stripping");
      vectorOfPreprocessingActionsAndNames[i].action->setDisabled(true);
      connect(vectorOfPreprocessingActionsAndNames[i].action, SIGNAL(triggered()), this, SLOT(ImageSkullStripping()));
    }
    //else if (vectorOfPreprocessingActionsAndNames[i].name.find("DCM2NIfTI") != std::string::npos)
    //{
    //  vectorOfPreprocessingActionsAndNames[i].action->setText("DICOM to NIfTI");
    //  connect(vectorOfPreprocessingActionsAndNames[i].action, SIGNAL(triggered()), this, SLOT(DCM2NIfTIConversion()));
    //}
    else if (vectorOfPreprocessingActionsAndNames[i].name.find("deepmedic_brain") != std::string::npos)
    {
      vectorOfPreprocessingActionsAndNames[i].action->setText("Skull Stripping (DeepLearning)"); // TBD set at source
      connect(vectorOfPreprocessingActionsAndNames[i].action, &QAction::triggered, this, [this] { ApplicationDeepMedicSegmentation(fDeepMedicDialog::SkullStripping); });
    }
    else if (vectorOfPreprocessingActionsAndNames[i].name.find("breastNormalize") != std::string::npos)
    {
      vectorOfPreprocessingActionsAndNames[i].action->setText("Mammogram Preprocessing");
      connect(vectorOfPreprocessingActionsAndNames[i].action, SIGNAL(triggered()), this, SLOT(ImageMamogramPreprocess()));
    }
    else if (vectorOfPreprocessingActionsAndNames[i].name.find("BraTSPipeline") != std::string::npos)
    {
      vectorOfPreprocessingActionsAndNames[i].action->setText("BraTS Pipeline");
      connect(vectorOfPreprocessingActionsAndNames[i].action, SIGNAL(triggered()), this, SLOT(ImageBraTSPipeline()));
    }
  }

  // add a single function for all preprocessing steps, this function will check for the specific names and then initiate that algorithm
  for (size_t i = 0; i < vectorOfDeepLearningActionsAndNames.size(); i++)
  {
    if (vectorOfDeepLearningActionsAndNames[i].name.find("deepmedic_tumor") != std::string::npos)
    {
      vectorOfDeepLearningActionsAndNames[i].action->setText("Brain Tumor Segmentation"); // TBD set at source
      connect(vectorOfDeepLearningActionsAndNames[i].action, &QAction::triggered, this, [this] { ApplicationDeepMedicSegmentation(fDeepMedicDialog::Tumor); });
    }
    else if (vectorOfDeepLearningActionsAndNames[i].name.find("deepmedic_brain") != std::string::npos)
    {
      vectorOfDeepLearningActionsAndNames[i].action->setText("Skull Stripping"); // TBD set at source
      connect(vectorOfDeepLearningActionsAndNames[i].action, &QAction::triggered, this, [this] { ApplicationDeepMedicSegmentation(fDeepMedicDialog::SkullStripping); });
    }
  }

  //connect(&fetalbrainpanel, SIGNAL(skullstripfun()), this, SLOT(FetalBrain_SkullStripfunc()));
  //connect(&fetalbrainpanel, SIGNAL(drawlinear()), this, SLOT(FetalBrain_Predict()));
  //connect(&fetalbrainpanel, SIGNAL(TrainNewFetalModel(std::string, std::string)), this, SLOT(FetalBrain_TrainNewModel(const std::string &, const std::string &)));


  connect(m_imagesTable, SIGNAL(itemSelectionChanged()), this, SLOT(DisplayChanged()));
  connect(m_imagesTable, SIGNAL(itemClicked(QTableWidgetItem*)), this, SLOT(DisplayChanged(QTableWidgetItem*)));

  connect(imagesPanel, SIGNAL(sigImageTableSelectionChanged()), this, SLOT(DisplayChanged()));

  connect(windowSpinBox, SIGNAL(editingFinished()), this, SLOT(WindowLevelEdited()));
  connect(levelSpinBox, SIGNAL(editingFinished()), this, SLOT(WindowLevelEdited()));

  connect(presetComboBox, SIGNAL(currentIndexChanged(int)), this, SLOT(UpdateWindowLevel()));
  connect(thresholdSpinBox, SIGNAL(valueChanged(double)), this, SLOT(thresholdSpinBoxChanged()));
  thresholdSpinBox->setValue(20);



  //connect(tumorPanel, SIGNAL(UpdateRenderWindows()), this, SLOT(UpdateRenderWindows()));
  //connect(tumorPanel, SIGNAL(SetActiveLandmarksTypeSignal(int, int, int)), this, SLOT(SetActiveLandmarksType(int, int, int)));
  //connect(tumorPanel, SIGNAL(MoveSlicerCursor(double, double, double)), this, SLOT(MoveSlicerCursor(double, double, double)));
  //connect(tumorPanel, SIGNAL(helpClicked_Interaction(std::string)), this, SLOT(help_contextual(std::string)));

  connect(segmentationPanel, SIGNAL(helpClicked_SegmentationUsage(std::string)), this, SLOT(help_contextual(std::string)));
  connect(trainingPanel, SIGNAL(helpClicked_SegmentationUsage(std::string)), this, SLOT(help_contextual(std::string)));

  connect(drawingPanel, SIGNAL(clearMask(int)), this, SLOT(clearMask(int)));
  connect(drawingPanel, SIGNAL(CurrentBrushSizeChanged(int)), this, SLOT(ChangeBrushSize(int)));
  connect(drawingPanel, SIGNAL(UndoButtonClicked()), this, SLOT(UndoFunctionality()));
  //connect(drawingPanel, SIGNAL(FillButtonClicked(int)), this, SLOT(FillLabel(int)));
  connect(drawingPanel, SIGNAL(shapesButtonClicked(int)), this, SLOT(updateDrawMode(int)));
  connect(drawingPanel, SIGNAL(CurrentDrawingLabelChanged(int)), this, SLOT(updateDrawMode()));
  connect(drawingPanel, SIGNAL(CurrentMaskOpacityChanged(int)), this, SLOT(ChangeMaskOpacity(int)));
  connect(drawingPanel, SIGNAL(helpClicked_Interaction(std::string)), this, SLOT(help_contextual(std::string)));
  connect(drawingPanel, SIGNAL(sig_ChangeLabelValuesClicked(const std::string, const std::string)), this, SLOT(CallLabelValuesChange(const std::string, const std::string)));


  //connect(&recurrencePanel, SIGNAL(SubjectBasedRecurrenceEstimate(std::string, bool, bool, bool, bool)), this, SLOT(StartRecurrenceEstimate(const std::string &, bool, bool, bool, bool)));
  //connect(&recurrencePanel, SIGNAL(SubjectBasedExistingRecurrenceEstimate(std::string, std::string, bool, bool, bool, bool)), this, SLOT(LoadedSubjectExistingRecurrenceEstimate(const std::string &, const std::string &, bool, bool, bool, bool)));
  //connect(&recurrencePanel, SIGNAL(ExistingModelBasedRecurrenceEstimate(std::string, std::string, std::string, bool, bool, bool, bool)), this, SLOT(RecurrenceEstimateOnExistingModel(const std::string &, const std::string &, const std::string &, bool, bool, bool, bool)));
  //connect(&recurrencePanel, SIGNAL(TrainNewModel(std::string, std::string, bool, bool, bool, bool)), this, SLOT(TrainNewModelOnGivenData(const std::string &, const std::string &, bool, bool, bool, bool)));

  //connect(&pseudoPanel, SIGNAL(ExistingModelBasedPseudoprogressionEstimate(std::string, std::string, std::string, bool, bool, bool, bool)), this, SLOT(PseudoprogressionEstimateOnExistingModel(const std::string &, const std::string &, const std::string &, bool, bool, bool, bool)));
  //connect(&pseudoPanel, SIGNAL(TrainNewPseudoModel(std::string, std::string, bool, bool, bool, bool)), this, SLOT(TrainNewPseudoprogressionModelOnGivenData(const std::string &, const std::string &, bool, bool, bool, bool)));


  //connect(&survivalPanel, SIGNAL(SurvivalPredictionOnExistingModel(const std::string, const std::string, const std::string)), this, SLOT(CallForSurvivalPredictionOnExistingModelFromMain(const std::string, const std::string, const std::string)));
  //connect(&survivalPanel, SIGNAL(PrepareNewSurvivalPredictionModel(const std::string, const std::string)), this, SLOT(CallForNewSurvivalPredictionModelFromMain(const std::string, const std::string)));

  //connect(&egfrv3Panel, SIGNAL(EGFRvIIIPredictionOnExistingModel(const std::string, const std::string, const std::string)), this, SLOT(CallForEGFRvIIIPredictionOnExistingModelFromMain(const std::string, const std::string, const std::string)));
  //connect(&egfrv3Panel, SIGNAL(PrepareNewEGFRvIIIPredictionModel(const std::string, const std::string)), this, SLOT(CallForNewEGFRvIIIPredictionModelFromMain(const std::string, const std::string)));


  //connect(&msubtypePanel, SIGNAL(MolecularSubtypePredictionOnExistingModel(const std::string, const std::string, const std::string)), this, SLOT(CallForMolecularSubtypePredictionOnExistingModelFromMain(const std::string, const std::string, const std::string)));
  //connect(&msubtypePanel, SIGNAL(PrepareNewMolecularSubtypePredictionModel(const std::string, const std::string)), this, SLOT(CallForNewMolecularSubtypePredictionModelFromMain(const std::string, const std::string)));


  //connect(&skullStrippingPanel, SIGNAL(RunSkullStripping(const std::string, const std::string, const std::string, const std::string)), this, SLOT(CallImageSkullStripping(const std::string, const std::string, const std::string, const std::string)));
  //connect(&dcmConverter, SIGNAL(RunDICOMConverter(const std::string, const std::string)), this, SLOT(CallDCM2NIfTIConversion(const std::string, const std::string)));
  connect(&histoMatchPanel, SIGNAL(RunHistogramMatching(const std::string, const std::string, const std::string)), this, SLOT(CallImageHistogramMatching(const std::string, const std::string, const std::string)));
  //connect(&deepMedicNormPanel, SIGNAL(RunDeepMedicNormalizer(const std::string, const std::string, const std::string, const std::string, const std::string, const std::string, const std::string, bool)), this, SLOT(CallImageDeepMedicNormalizer(const std::string, const std::string, const std::string, const std::string, const std::string, const std::string, const std::string, bool)));
  //connect(&directionalityEstimator, SIGNAL(RunDirectionalityEstimator(const std::string, const std::string, const std::string)), this, SLOT(CallDirectionalityEstimator(const std::string, const std::string, const std::string)));
  connect(&bratsPipelineDialog, SIGNAL(RunBraTSPipeline(const std::string, const std::string, const std::string, const std::string, const std::string)), this, SLOT(CallBraTSPipeline(const std::string, const std::string, const std::string, const std::string, const std::string)));

  //connect(&pcaPanel, SIGNAL(ExistingModelBasedPCAEstimate(std::string, std::string, std::string)), this, SLOT(PCAEstimateOnExistingModel(const std::string &, const std::string &, const std::string &)));
  //connect(&pcaPanel, SIGNAL(TrainNewPCAModel(std::string, std::string)), this, SLOT(TrainNewPCAModelOnGivenData(const std::string &, const std::string &)));


  //connect(&pcaPanel, SIGNAL(RunPCAEstimation(const int, const std::string, const std::string)), this, SLOT(CallPCACalculation(const int, const std::string, const std::string)));
  //connect(&trainingPanel, SIGNAL(RunTrainingSimulation(const std::string, const std::string, const std::string, const std::string, int, int, int)), this, SLOT(CallTrainingSimulation(const std::string, const std::string, const std::string, const std::string, int, int, int)));

  //connect(&perfmeasuresPanel, SIGNAL(RunPerfusionMeasuresCalculation(const double, const bool, const bool, const bool, const std::string, const std::string)), this, SLOT(CallPerfusionMeasuresCalculation(const double, const bool, const bool, const bool, const std::string, const std::string)));
  //connect(&perfalignPanel, SIGNAL(RunPerfusionAlignmentCalculation(double,int, int,const std::string, const std::string, const std::string, const std::string)), this, SLOT(CallPerfusionAlignmentCalculation(double,int, int, const std::string, const std::string, const std::string, const std::string)));


  //connect(&diffmeasuresPanel, SIGNAL(RunDiffusionMeasuresCalculation(const std::string, const std::string, const std::string, const std::string, const bool, const bool, const bool, const bool, const std::string)), this,
    //SLOT(CallDiffusionMeasuresCalculation(const std::string, const std::string, const std::string, const std::string, const bool, const bool, const bool, const bool, const std::string)));

  //connect(&whiteStripeNormalizer, SIGNAL(RunWhiteStripe(double, int, int, int, double, double, int, bool, const std::string)), this, SLOT(CallWhiteStripe(double, int, int, int, double, double, int, bool, const std::string)));

  //connect(&atlasPanel, SIGNAL(GeneratePopualtionAtlas(const std::string, const std::string, const std::string, const std::string)), this, SLOT(CallGeneratePopualtionAtlas(const std::string, const std::string, const std::string, const std::string)));

  //connect(&nodulePanel, SIGNAL(SBRTNoduleParamReady(const std::string, const int)), this, SLOT(CallSBRTNodule(const std::string, const int)));

  connect(&deepMedicDialog, SIGNAL(RunDeepMedic(const std::string, const std::string)), this, SLOT(CallDeepMedicSegmentation(const std::string, const std::string)));
  //connect(&texturePipelineDialog, SIGNAL(RunTextureFeaturePipeline(const std::string)), this, SLOT(CallTexturePipeline(const std::string)));

  //connect(this, SIGNAL(SeedPointsFocused(bool)), tumorPanel, SLOT(sTableFocused(bool)));
  //connect(this, SIGNAL(TissuePointsFocused(bool)), tumorPanel, SLOT(tTableFocused(bool)));
  connect(m_tabWidget, SIGNAL(currentChanged(int)), this, SLOT(panelChanged(int)));
  connect(infoPanel, SIGNAL(MoveSlicerCursor(double, double, double, int)), this, SLOT(MoveSlicerCursor(double, double, double, int)));

  AxialViewWidget->hide();
  CoronalViewWidget->hide();
  SaggitalViewWidget->hide();
  infoPanel->show();

  windowLabel->setEnabled(false);
  windowSpinBox->setEnabled(false);
  levelSpinBox->setEnabled(false);
  levelLabel->setEnabled(false);

  thresholdLabel->setEnabled(false);
  thresholdSpinBox->setEnabled(false);
  presetLabel->setEnabled(false);
  presetComboBox->setEnabled(false);

  setAcceptDrops(true);

  mCustomImageToThreshold_min = 0;
  mCustomImageToThreshold_max = 0;

  m_drawShapeMode = SHAPE_MODE_NONE;

  m_messageLabel = new QLabel(":");
  statusBar()->addPermanentWidget(m_messageLabel);
  m_progressBar = new QProgressBar(this);
  QSize tempSize = this->size();
  m_progressBar->setFixedWidth(tempSize.width() * 0.75);
  statusBar()->addPermanentWidget(m_progressBar);
  m_progressBar->setValue(0);

  mHelpDlg = new fHelpDialog();

  //connect
  connect(m_toolTabdock, SIGNAL(topLevelChanged(bool)), this, SLOT(toolTabDockChanged(bool)));

  //recurrencePanel.SetCurrentLoggerPath(m_tempFolderLocation);
  //msubtypePanel.SetCurrentLoggerPath(m_tempFolderLocation);
  //survivalPanel.SetCurrentLoggerPath(m_tempFolderLocation);

  //
  actionLoad_Nifti_Images->setText(QApplication::translate("fMainWindow", "Image(s)", 0));
  actionLoad_Nifti_ROI->setText(QApplication::translate("fMainWindow", "ROI", 0));
  actionLoad_Dicom_Images->setText(QApplication::translate("fMainWindow", "Dicom", 0));

  actionPreferences->setText(QApplication::translate("fMainWindow", "Preferences", 0));
  actionSave_Nifti_Images->setText(QApplication::translate("fMainWindow", "Image (NIfTI)", 0));
  actionSave_Dicom_Images->setText(QApplication::translate("fMainWindow", "Image (DICOM)", 0));
  actionSave_ROI_Images->setText(QApplication::translate("fMainWindow", "ROI (NIfTI)", 0));
  actionSave_ROI_Dicom_Images->setText(QApplication::translate("fMainWindow", "ROI (DICOM)", 0));
  actionHelp_Interactions->setText(QApplication::translate("fMainWindow", "Usage", 0));
  help_discussion->setText(QApplication::translate("fMainWindow", "Discussion Forum", 0));
  help_forum->setText(QApplication::translate("fMainWindow", "Help Forum", 0));
  help_bugs->setText(QApplication::translate("fMainWindow", "Bugs and Feature", 0));
  help_features->setText(QApplication::translate("fMainWindow", "Feature Requests", 0));
  help_download->setText(QApplication::translate("fMainWindow", "Latest Downloads", 0));
  actionAbout->setText(QApplication::translate("fMainWindow", "About", 0));
  actionExit->setText(QApplication::translate("fMainWindow", "Exit", 0));
  actionAppGeodesic->setText(QApplication::translate("fMainWindow", "Geodesic segmentation", 0));
  actionAppGeodesicTraining->setText(QApplication::translate("fMainWindow", "Geodesic Training Segmentation", 0));
  //m_tabWidget->setTabText(m_tabWidget->indexOf(tumorPanel), QApplication::translate("fMainWindow", "Seed Points", 0));
  m_tabWidget->setTabText(m_tabWidget->indexOf(drawingPanel), QApplication::translate("fMainWindow", "Drawing", 0));
  m_tabWidget->setTabText(m_tabWidget->indexOf(imagesPanel), QApplication::translate("fMainWindow", "Images", 0));
  m_tabWidget->setTabText(m_tabWidget->indexOf(segmentationPanel), QApplication::translate("fMainWindow", "Segmentation", 0));
  m_tabWidget->setTabText(m_tabWidget->indexOf(trainingPanel), QApplication::translate("fMainWindow", "Training", 0));

  connect(segmentationPanel, SIGNAL(m_btnComputeClicked()), this, SLOT(OnSegmentationClicked()));
  connect(trainingPanel, SIGNAL(m_btnComputeClicked()), this, SLOT(OnTrainingClicked()));
}

fMainWindow::~fMainWindow()
{
  for (int i = 0; i < (int)mSlicerManagers.size(); i++)
  {
    if (mSlicerManagers[i] != NULL)
    {
      delete mSlicerManagers[i];
    }
  }
  if (mLandmarks)
  {
    delete mLandmarks;
  }
  if (mSeedPoints)
  {
    delete mSeedPoints;
  }
  if (mTissuePoints)
  {
    delete mTissuePoints;
  }

  // delete the temp directory every single time
  //if (cbica::filesInDirectory(m_tempFolderLocation).empty())
  {
    cbica::deleteDir(m_tempFolderLocation);
  }

  if (m_skipTutorialOnNextRun)
  {
    std::ofstream file;
    file.open(tutorialScreen.c_str());
    file << "User doesn't need the tutorial screen.\n";
    file.close();
  }

  if (mHelpDlg)
    delete mHelpDlg;

  if (segmentationPanel)
  {
    delete segmentationPanel;
  }

  if (trainingPanel)
  {
    delete trainingPanel;
  }

  ApplicationPreferences::GetInstance()->SerializePreferences();
}

void fMainWindow::ImageBraTSPipeline()
{
  // open a simple dialog box with reference image, input and output
  bratsPipelineDialog.SetCurrentImagePath(mInputPathName);
  bratsPipelineDialog.exec();
}

void fMainWindow::CallBraTSPipeline(const std::string t1ceImage, const std::string t1Image, const std::string t2Image, const std::string flImage, const std::string outputDir)
{
  if (!t1ceImage.empty() && !t1Image.empty() && !t2Image.empty() && !flImage.empty() && !outputDir.empty())
  {
    auto bratsPipelineExe = getApplicationPath("BraTSPipeline");
    if (!cbica::exists(bratsPipelineExe))
    {
      ShowErrorMessage("Could not find the BraTSPipeline executable");
      return;
    }

    QStringList args;
    args << "-t1" << t1Image.c_str() <<
      "-t1c" << t1ceImage.c_str() <<
      "-t2" << t2Image.c_str() <<
      "-fl" << flImage.c_str() <<
      "-o" << outputDir.c_str();

    QMessageBox *box = new QMessageBox(QMessageBox::Question, "Long running Application",
      "BraTS Pipeline takes ~30 minutes to run, during which FeTS UI will not be responsive; press OK to continue...",
      QMessageBox::Ok | QMessageBox::Cancel);
    box->setAttribute(Qt::WA_DeleteOnClose); //makes sure the msgbox is deleted automatically when closed
    box->setWindowModality(Qt::NonModal);
    QCoreApplication::processEvents();
    if (box->exec() == QMessageBox::Ok)
    {
      updateProgress(5, "Starting BraTS Pipeline");

      if (startExternalProcess(bratsPipelineExe.c_str(), args) != 0)
      {
        ShowErrorMessage("BraTS Pipeline returned with exit code != 0");
        updateProgress(0, "");
        return;
      }
    }
  }
  else
  {
    ShowErrorMessage("All input images need to be provided for BraTS Pipeline to run.");
    return;
  }
}

  void fMainWindow::loadFromCommandLine(std::vector< QString > files, bool comparisonMode, const std::string &maskImage, const float maskOpacity,
    const std::string &tumorPointFile, const std::string &tissuePointFile, bool firstRun)
  {
    auto qvectorString = QVector< QString >::fromStdVector(files);
    auto lst = QStringList::fromVector(QVector< QString >::fromStdVector(files));
    this->openImages(lst, true);
    if (!maskImage.empty())
    {
      this->readMaskFile(maskImage);
      this->ChangeMaskOpacity(maskOpacity * 10);
    }
    if (!tumorPointFile.empty())
    {
      //this->tumorPanel->sLoad(tumorPointFile.c_str());
    }
    if (!tissuePointFile.empty())
    {
      //this->tumorPanel->tLoad(tissuePointFile.c_str());
    }
    if (comparisonMode)
    {
      this->imagesPanel->CompareButtonClick();
    }

#ifdef CAPTK_PACKAGE_PROJECT
    if (firstRun)
    {
      this->CloseAllImages();
    }
#endif
  }

std::string fMainWindow::ConversionFrom2Dto3D(const std::string &fileName)
{
  using ImageTypeFloat2D = itk::Image< float, 2 >;
  auto reader = itk::ImageFileReader< ImageTypeFloat2D >::New();
  reader->SetFileName(fileName);
  auto ext = cbica::getFilenameExtension(fileName);
  if (cbica::IsDicom(fileName))
  {
    reader->SetImageIO(itk::GDCMImageIO::New());
  }
  else if ((ext == ".nii") || (ext == ".nii.gz"))
  {
    reader->SetImageIO(itk::NiftiImageIO::New());
  }

  try
  {
    reader->Update();
  }
  catch (itk::ExceptionObject& e)
  {
    ShowErrorMessage("Exception caught while reading the image '" + fileName + "':\n\n" + e.what());
    return "";
  }

  auto image_2D = reader->GetOutput();
  auto index2D = image_2D->GetLargestPossibleRegion().GetIndex();
  auto size2D = image_2D->GetLargestPossibleRegion().GetSize();
  auto spacing2D = image_2D->GetSpacing();
  auto origin2D = image_2D->GetOrigin();

  // write into m_tempFolderLocation and then read pass that file to LoadSlicerImages
  auto image_3D = ImageTypeFloat3D::New();
  ImageTypeFloat3D::RegionType region;
  ImageTypeFloat3D::RegionType::SizeType size;
  ImageTypeFloat3D::RegionType::IndexType start;
  ImageTypeFloat3D::PointType origin;
  ImageTypeFloat3D::SpacingType spacing;

  // populate the region to place the slice in
  size[0] = size2D[0];
  size[1] = size2D[1];
  size[2] = 1;
  start[0] = index2D[0];
  start[1] = index2D[1];
  start[2] = 0;
  spacing[0] = image_2D->GetSpacing()[0];
  spacing[1] = image_2D->GetSpacing()[1];
  spacing[2] = 1;
  origin[0] = origin2D[0];
  origin[1] = origin2D[1];
  origin[2] = 0;

  region.SetSize(size);
  region.SetIndex(start);
  image_3D->SetRegions(region);
  image_3D->SetRequestedRegion(region);
  image_3D->SetBufferedRegion(region);
  image_3D->Allocate();
  image_3D->FillBuffer(0);
  image_3D->SetOrigin(origin);
  image_3D->SetSpacing(spacing);

  itk::ImageRegionIteratorWithIndex <ImageTypeFloat3D> iter(image_3D, image_3D->GetLargestPossibleRegion());
  itk::ImageRegionIteratorWithIndex <ImageTypeFloat2D> iter2d(image_2D, image_2D->GetLargestPossibleRegion());

  for (iter.GoToBegin(), iter2d.GoToBegin(); !iter2d.IsAtEnd(); ++iter, ++iter2d)
  {
    iter.Set(iter2d.Get());
  }

  auto imageName = m_tempFolderLocation + "/" + cbica::getFilenameBase(fileName) + ".nii.gz";
  cbica::WriteImage< ImageTypeFloat3D >(image_3D, imageName);

  return imageName;
}

void fMainWindow::about()
{
#if CAPTK_PACKAGE_PROJECT
  mHelpTutorial.exec();
#endif
}

void fMainWindow::help_Interactions()
{
  openLink("https://fets-ai.github.io/Front-End/");
}

void fMainWindow::help_Download(QAction* action)
{
  auto currentApp = action->text().toStdString();
  std::string path = getCaPTkDataDir();
  auto currentLink = "ftp://www.nitrc.org/home/groups/captk/downloads/SampleData_1.6.0/" + currentApp + ".zip";
  cbica::Logging(loggerFile, currentLink);
  if (!openLink(currentLink))
  {
      ShowErrorMessage("CaPTk couldn't open the browser to download specified sample data.", this);
    return;
  }
}

void fMainWindow::help_BugTracker()
{
  if (!openLink("https://github.com/CBICA/CaPTk/issues"))
  {
    ShowErrorMessage("CaPTk couldn't open the browser to open the Bug Tracker");
    return;
  }
}

void fMainWindow::EnableThresholdOfMask()
{

  // only do calculations on current image(s)
  auto items = m_imagesTable->selectedItems();
  if (items.empty())
  {
    return;
  }
  int index = GetSlicerIndexFromItem(items[0]);
  if (index < 0 || index >= (int)mSlicerManagers.size())
  {
    return;
  }

  typedef itk::MinimumMaximumImageCalculator < ImageTypeFloat3D > ImageCalculatorFilterType;
  ImageCalculatorFilterType::Pointer imageCalculatorFilter = ImageCalculatorFilterType::New();
  imageCalculatorFilter->SetImage(mSlicerManagers[index]->mITKImage);
  imageCalculatorFilter->Compute();
  double actualMin = imageCalculatorFilter->GetMinimum();
  double actualMax = imageCalculatorFilter->GetMaximum();

  thresholdSpinBox->setRange(actualMin, actualMax);
  thresholdSpinBox->setValue((actualMin + actualMax) / 2);
}

void fMainWindow::SaveImage_withFile(int indexOfInputImageToWrite, QString saveFileName)
{
  auto index = indexOfInputImageToWrite;
  if (!saveFileName.isEmpty())
  {
    auto saveFileName_string = saveFileName.toStdString();
    typedef ImageTypeFloat3D ImageType;
    ImageType::DirectionType originalDirection;
    originalDirection[0][0] = mSlicerManagers[index]->mDirection(0, 0);
    originalDirection[0][1] = mSlicerManagers[index]->mDirection(0, 1);
    originalDirection[0][2] = mSlicerManagers[index]->mDirection(0, 2);
    originalDirection[1][0] = mSlicerManagers[index]->mDirection(1, 0);
    originalDirection[1][1] = mSlicerManagers[index]->mDirection(1, 1);
    originalDirection[1][2] = mSlicerManagers[index]->mDirection(1, 2);
    originalDirection[2][0] = mSlicerManagers[index]->mDirection(2, 0);
    originalDirection[2][1] = mSlicerManagers[index]->mDirection(2, 1);
    originalDirection[2][2] = mSlicerManagers[index]->mDirection(2, 2);

    ImageType::PointType originalOrigin;
    originalOrigin = mSlicerManagers[index]->mOrigin;

    if (mSlicerManagers[index]->GetPreset() == PRESET_THRESHOLD)
    {
      auto img = convertVtkToItk<ImageTypeFloat3D::PixelType, ImageTypeFloat3D::ImageDimension>(mSlicerManagers[index]->mImage);
      double threshold = mSlicerManagers[index]->mThreshold;
      //
      auto duplicator = itk::ImageDuplicator< ImageTypeFloat3D >::New();
      duplicator->SetInputImage(img);
      duplicator->Update();
      auto seg = duplicator->GetOutput();
      //
      itk::ImageRegionIterator< ImageTypeFloat3D > imgIterator(seg, seg->GetLargestPossibleRegion());
      for (imgIterator.GoToBegin(); !imgIterator.IsAtEnd(); ++imgIterator)
      {
        // TBD: this was the original action to handle threshold -- needs testing
        if (imgIterator.Get() <= threshold)
        {
          imgIterator.Set(1);
        }
        else
        {
          imgIterator.Set(0);
        }
      }
      //

      auto infoChanger = itk::ChangeInformationImageFilter< ImageType >::New();
      infoChanger->SetInput(seg);
      infoChanger->ChangeDirectionOn();
      infoChanger->ChangeOriginOn();
      infoChanger->SetOutputDirection(originalDirection);
      infoChanger->SetOutputOrigin(originalOrigin);
      infoChanger->Update();

      cbica::WriteImage< ImageTypeFloat3D >(infoChanger->GetOutput(), correctExtension(saveFileName_string));
    }
    else
    {
      auto img = convertVtkToItk< ImageType::PixelType, ImageTypeFloat3D::ImageDimension>(mSlicerManagers[index]->mImage);

      auto infoChanger = itk::ChangeInformationImageFilter< ImageType >::New();
      infoChanger->SetInput(img);
      infoChanger->ChangeDirectionOn();
      infoChanger->ChangeOriginOn();
      infoChanger->SetOutputDirection(originalDirection);
      infoChanger->SetOutputOrigin(originalOrigin);
      infoChanger->Update();

      cbica::WriteImage< ImageType >(infoChanger->GetOutput(), correctExtension(saveFileName_string));

      std::string InputPixelType = mSlicerManagers[index]->mImage->GetScalarTypeAsString();
      if (InputPixelType == "short")
      {
        using ImageTypeToWrite = itk::Image<short, ImageTypeFloat3D::ImageDimension>;
        auto img = convertVtkToItk<ImageTypeToWrite::PixelType, ImageTypeFloat3D::ImageDimension>(mSlicerManagers[index]->mImage);
        auto infoChanger = itk::ChangeInformationImageFilter< ImageTypeToWrite >::New();
        infoChanger->SetInput(img);
        infoChanger->ChangeDirectionOn();
        infoChanger->ChangeOriginOn();
        infoChanger->SetOutputDirection(originalDirection);
        infoChanger->SetOutputOrigin(originalOrigin);
        infoChanger->Update();

        cbica::WriteImage< ImageTypeToWrite >(infoChanger->GetOutput(), correctExtension(saveFileName_string));
      }
      else if (InputPixelType == "unsigned short")
      {
        using ImageTypeToWrite = itk::Image<unsigned short, ImageTypeFloat3D::ImageDimension>;
        auto img = convertVtkToItk<ImageTypeToWrite::PixelType, ImageTypeFloat3D::ImageDimension>(mSlicerManagers[index]->mImage);
        auto infoChanger = itk::ChangeInformationImageFilter< ImageTypeToWrite >::New();
        infoChanger->SetInput(img);
        infoChanger->ChangeDirectionOn();
        infoChanger->ChangeOriginOn();
        infoChanger->SetOutputDirection(originalDirection);
        infoChanger->SetOutputOrigin(originalOrigin);
        infoChanger->Update();

        cbica::WriteImage< ImageTypeToWrite >(infoChanger->GetOutput(), correctExtension(saveFileName_string));
      }
      else if (InputPixelType == "char")
      {
        using ImageTypeToWrite = itk::Image<char, ImageTypeFloat3D::ImageDimension>;
        auto img = convertVtkToItk<ImageTypeToWrite::PixelType, ImageTypeFloat3D::ImageDimension>(mSlicerManagers[index]->mImage);
        auto infoChanger = itk::ChangeInformationImageFilter< ImageTypeToWrite >::New();
        infoChanger->SetInput(img);
        infoChanger->ChangeDirectionOn();
        infoChanger->ChangeOriginOn();
        infoChanger->SetOutputDirection(originalDirection);
        infoChanger->SetOutputOrigin(originalOrigin);
        infoChanger->Update();

        cbica::WriteImage< ImageTypeToWrite >(infoChanger->GetOutput(), correctExtension(saveFileName_string));
      }
      else if (InputPixelType == "unsigned char")
      {
        using ImageTypeToWrite = itk::Image<unsigned char, ImageTypeFloat3D::ImageDimension>;
        auto img = convertVtkToItk<ImageTypeToWrite::PixelType, ImageTypeFloat3D::ImageDimension>(mSlicerManagers[index]->mImage);
        auto infoChanger = itk::ChangeInformationImageFilter< ImageTypeToWrite >::New();
        infoChanger->SetInput(img);
        infoChanger->ChangeDirectionOn();
        infoChanger->ChangeOriginOn();
        infoChanger->SetOutputDirection(originalDirection);
        infoChanger->SetOutputOrigin(originalOrigin);
        infoChanger->Update();

        cbica::WriteImage< ImageTypeToWrite >(infoChanger->GetOutput(), correctExtension(saveFileName_string));
      }
      else if (InputPixelType == "int")
      {
        using ImageTypeToWrite = itk::Image<int, ImageTypeFloat3D::ImageDimension>;
        auto img = convertVtkToItk<ImageTypeToWrite::PixelType, ImageTypeFloat3D::ImageDimension>(mSlicerManagers[index]->mImage);
        auto infoChanger = itk::ChangeInformationImageFilter< ImageTypeToWrite >::New();
        infoChanger->SetInput(img);
        infoChanger->ChangeDirectionOn();
        infoChanger->ChangeOriginOn();
        infoChanger->SetOutputDirection(originalDirection);
        infoChanger->SetOutputOrigin(originalOrigin);
        infoChanger->Update();

        cbica::WriteImage< ImageTypeToWrite >(infoChanger->GetOutput(), correctExtension(saveFileName_string));
      }
      else if (InputPixelType == "unsigned int")
      {
        using ImageTypeToWrite = itk::Image<unsigned int, ImageTypeFloat3D::ImageDimension>;
        auto img = convertVtkToItk<ImageTypeToWrite::PixelType, ImageTypeFloat3D::ImageDimension>(mSlicerManagers[index]->mImage);
        auto infoChanger = itk::ChangeInformationImageFilter< ImageTypeToWrite >::New();
        infoChanger->SetInput(img);
        infoChanger->ChangeDirectionOn();
        infoChanger->ChangeOriginOn();
        infoChanger->SetOutputDirection(originalDirection);
        infoChanger->SetOutputOrigin(originalOrigin);
        infoChanger->Update();

        cbica::WriteImage< ImageTypeToWrite >(infoChanger->GetOutput(), correctExtension(saveFileName_string));
      }
      else if (InputPixelType == "double")
      {
        using ImageTypeToWrite = itk::Image<double, ImageTypeFloat3D::ImageDimension>;
        auto img = convertVtkToItk<ImageTypeToWrite::PixelType, ImageTypeFloat3D::ImageDimension>(mSlicerManagers[index]->mImage);
        auto infoChanger = itk::ChangeInformationImageFilter< ImageTypeToWrite >::New();
        infoChanger->SetInput(img);
        infoChanger->ChangeDirectionOn();
        infoChanger->ChangeOriginOn();
        infoChanger->SetOutputDirection(originalDirection);
        infoChanger->SetOutputOrigin(originalOrigin);
        infoChanger->Update();

        cbica::WriteImage< ImageTypeToWrite >(infoChanger->GetOutput(), correctExtension(saveFileName_string));
      }
      else if (InputPixelType == "float")
      {
        using ImageTypeToWrite = itk::Image<float, ImageTypeFloat3D::ImageDimension>;
        auto img = convertVtkToItk<ImageTypeToWrite::PixelType, ImageTypeFloat3D::ImageDimension>(mSlicerManagers[index]->mImage);
        auto infoChanger = itk::ChangeInformationImageFilter< ImageTypeToWrite >::New();
        infoChanger->SetInput(img);
        infoChanger->ChangeDirectionOn();
        infoChanger->ChangeOriginOn();
        infoChanger->SetOutputDirection(originalDirection);
        infoChanger->SetOutputOrigin(originalOrigin);
        infoChanger->Update();

        cbica::WriteImage< ImageTypeToWrite >(infoChanger->GetOutput(), correctExtension(saveFileName_string));
      }
      else
      {
        cbica::Logging(loggerFile, "Error, input pixel type, '" + InputPixelType + "' is unknown!");
      }
      updateProgress(0, "Image saved! (" + saveFileName_string + ")");
    }
  }
}

void fMainWindow::SaveImage()
{
  auto items = m_imagesTable->selectedItems();
  if (items.empty()) 
  {
    return;
  }
  int index = GetSlicerIndexFromItem(items[0]);
  if (index < 0 || index >= (int)mSlicerManagers.size()) 
  {
    return;
  }
  //
  QString saveFileName = getSaveFile(this, mInputPathName, mInputPathName + "_new.nii.gz");
  SaveImage_withFile(index, saveFileName);
}

void fMainWindow::InitMask(vtkImageData* image)
{
  int extent[6];
  double spacing[3];
  double origin[3];
  int i, j, k;

  image->GetExtent(extent);
  image->GetSpacing(spacing);
  image->GetOrigin(origin);

  mMask->Initialize();
  mMask->SetExtent(extent);
  mMask->SetSpacing(spacing);
  mMask->SetOrigin(origin);
#if VTK_MAJOR_VERSION <= 5
  mMask->SetScalarTypeToFloat();
  mMask->SetNumberOfScalarComponents(1);
  mMask->AllocateScalars();
#else
  mMask->AllocateScalars(VTK_FLOAT, 1);
#endif
  {
    int vd_x, vd_y, vd_z;
    vd_x = mMask->GetDimensions()[0];
    vd_y = mMask->GetDimensions()[1];
    vd_z = mMask->GetDimensions()[2];
    float* pData = (float*)mMask->GetScalarPointer();
    for (k = 0; k < vd_z; k++)
    {
      for (j = 0; j < vd_y; j++)
      {
        for (i = 0; i < vd_x; i++)
        {
          *pData = 0;
          pData++;
        }
      }
    }
  }
}
void fMainWindow::updateDrawMode(int shapeMode)
{
  if (shapeMode >= 0)
  {
    m_drawShapeMode = SHAPE_MODE(shapeMode);
  }
  if (m_drawShapeMode == SHAPE_MODE_NONE)
  {
    AxialViewWidget->unsetCursor();
    CoronalViewWidget->unsetCursor();
    SaggitalViewWidget->unsetCursor();
    return;
  }

  int color[4] = { 0, 0, 0, 0 };
  int drawLabel = this->getSelectedDrawLabel();
  int drawSize = this->getSelectedDrawSize();

  switch (drawLabel)
  {
  case DRAW_MODE_LABEL_1: // near
    color[0] = 255;
    color[1] = 0;
    color[2] = 0;
    color[3] = 255;
    break;
  case DRAW_MODE_LABEL_2: // far
    color[0] = 0;
    color[1] = 255;
    color[2] = 0;
    color[3] = 255;
    break;
  case DRAW_MODE_LABEL_3:
    color[0] = 255;
    color[1] = 255;
    color[2] = 0;
    color[3] = 255;
    break;
  case DRAW_MODE_LABEL_4:
    color[0] = 0;
    color[1] = 0;
    color[2] = 255;
    color[3] = 255;
    break;
  case DRAW_MODE_LABEL_5:
    color[0] = 255;
    color[1] = 0;
    color[2] = 255;
    color[3] = 255;
    break;


  default:
    break;
  }
  if (m_drawShapeMode == SHAPE_MODE_ERASER)
  {
    color[0] = 255;
    color[1] = 255;
    color[2] = 255;
    color[3] = 255;
  }

  QImage img(33, 33, QImage::Format_ARGB32);
  for (int j = -16; j <= 16; j++)
  {
    for (int i = -16; i <= 16; i++)
    {
      if (abs(i) <= drawSize && abs(j) <= drawSize)
      {
        img.setPixel(i + 16, j + 16, qRgba(color[0], color[1], color[2], color[3]));
      }
      else
      {
        img.setPixel(i + 16, j + 16, qRgba(0, 0, 0, 0));
      }
    }
  }
  AxialViewWidget->setCursor(QCursor(Qt::CrossCursor));
  CoronalViewWidget->setCursor(QCursor(Qt::CrossCursor));
  SaggitalViewWidget->setCursor(QCursor(Qt::CrossCursor));


}

void fMainWindow::LoadNonViewingImages(const std::string &directoryname, const int &imagetype_int, const int &imagesubtype)
{
  for (unsigned int index = 0; index < mNonViewingImageManager.size(); index++)
  {
    if (mNonViewingImageManager[index]->GetPathFileName() == directoryname)
    {
      return;
    }
  }
  SimpleImageManager* nonViewingImage = new SimpleImageManager();
  nonViewingImage->SetPathFileName(directoryname);
  nonViewingImage->SetFileName(directoryname);
  if (imagetype_int == CAPTK::NIfTI)
    nonViewingImage->ReadGivenNonViewingNiftiImage(directoryname, imagesubtype);
  else
    nonViewingImage->ReadGivenNonViewingDicomImage(directoryname, imagesubtype);

  if (nonViewingImage->GetLastEncounteredError() != "")
    delete nonViewingImage;
  else
    mNonViewingImageManager.push_back(nonViewingImage);


  //Add image Information in NonViewing images table
  //----------------------------------------------
  nonViewingImage->mImageType = imagetype_int;
  nonViewingImage->mImageSubType = imagesubtype;

  std::string strImageType;
  if (nonViewingImage->mImageSubType == CAPTK::ImageModalityType::IMAGE_TYPE_DTI)
    strImageType = "DTI";

  int rowindex = m_nonVisImagesTable->rowCount();
  m_nonVisImagesTable->setRowCount(rowindex + 1);

  QFileInfo fileinfo(nonViewingImage->GetFileName().c_str());
  QString id = directoryname.c_str() + QString::number(mNonViewingImageManager.size() - 1);
  {
    QTableWidgetItem *item = new QTableWidgetItem(fileinfo.fileName());
    item->setData(Qt::UserRole, directoryname.c_str());

    QTablePushButton* cButton = new QTablePushButton;
    cButton->setItem(item);
    cButton->setText(tr("X"));
    cButton->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);


    QString in = QString::fromStdString(strImageType);
    QTableWidgetItem *item2 = new QTableWidgetItem(in);
    item2->setData(Qt::UserRole, in.toStdString().c_str());

    if (nonViewingImage->mImageSubType == CAPTK::ImageModalityType::IMAGE_TYPE_DTI)
      connect(cButton, SIGNAL(clickedInto(QTableWidgetItem*)), this, SLOT(CloseNonViewingDTIImage(QTableWidgetItem*)));

    m_nonVisImagesTable->setCellWidget(rowindex, TAB_IMAGES_COLUMN_CLOSE, cButton);
    m_nonVisImagesTable->setItem(rowindex, TAB_IMAGES_COLUMN_NAME, item);
    m_nonVisImagesTable->setItem(rowindex, TAB_IMAGES_COLUMN_TYPE, item2);

    m_nonVisImagesTable->resizeRowsToContents();
  }
}



void fMainWindow::LoadSlicerImages(const std::string &fileName, const int &imagetype_int, bool bSkipDup)
{
  std::string fname;
  auto extension = cbica::getFilenameExtension(fileName);
  //if (extension != ".dcm")
  {
    if (extension == ".zip")
    {
      ShowErrorMessage("Please extract the zip file before trying to load into CaPTk.");
      return;
    }
    if ((extension != ".nii") && (extension != ".nii.gz"))
    {
      ShowErrorMessage("Only DICOM (dcm) or NIfTI (nii/nii.gz) images are supported right now; please contact CBICA for adding extended support");
      return;
    }
    //if ((extension == ".dcm") || 
    //  (extension == ".DCM") ||
    //  (extension == ".dicom") || 
    //  (extension == "") || 
    //  (extension == ".ima") ||
    //  (extension == ".IMA"))
    if (cbica::IsDicom(fileName))
    {
      QDir d = QFileInfo(fileName.c_str()).absoluteDir();
      fname = d.absolutePath().toStdString();
      dicomfilename = fname;
    }
    else 
      fname = fileName;
    auto imageInfo = cbica::ImageInfo(fname);
    SlicerManager* imageManager = new SlicerManager(3, mLandmarks, mSeedPoints, mTissuePoints);
    imageManager->mImageSubType = CAPTK::ImageModalityType::IMAGE_TYPE_UNDEFINED;
    imageManager->SetComparisonMode(false);

    bool bFirstLoad = false;
    if (mSlicerManagers.empty())
    {
      bFirstLoad = true;
    }
    if (imageInfo.GetImageDimensions() == 2)
    {
      fname = ConversionFrom2Dto3D(fname);
    }
    else if (!bFirstLoad)
    {
      {
        //auto temp_prev = cbica::normPath(m_tempFolderLocation + "/temp_prev.nii.gz");
        //ImageTypeFloat3D::DirectionType originaldirection;
        //originaldirection[0][0] = mSlicerManagers[0]->mDirection(0, 0);
        //originaldirection[0][1] = mSlicerManagers[0]->mDirection(0, 1);
        //originaldirection[0][2] = mSlicerManagers[0]->mDirection(0, 2);
        //originaldirection[1][0] = mSlicerManagers[0]->mDirection(1, 0);
        //originaldirection[1][1] = mSlicerManagers[0]->mDirection(1, 1);
        //originaldirection[1][2] = mSlicerManagers[0]->mDirection(1, 2);
        //originaldirection[2][0] = mSlicerManagers[0]->mDirection(2, 0);
        //originaldirection[2][1] = mSlicerManagers[0]->mDirection(2, 1);
        //originaldirection[2][2] = mSlicerManagers[0]->mDirection(2, 2);
        //auto img = convertVtkToItk< ImageTypeFloat3D::PixelType, ImageTypeFloat3D::ImageDimension >(mSlicerManagers[0]->mImage);
        //img->SetDirection(originaldirection);

        //cbica::WriteImage< ImageTypeFloat3D >(img, temp_prev);

        bool fourDImage = false;
        if ((imageInfo.GetImageDimensions() == 4) || (mSlicerManagers[0]->mImageSubType == CAPTK::ImageModalityType::IMAGE_TYPE_PERFUSION))
        {
          fourDImage = true;
        }
        if (!cbica::ImageSanityCheck(fname, mSlicerManagers[0]->GetPathFileName(), fourDImage))
        {
          ShowErrorMessage("The physical dimensions of the previously loaded image and current image are inconsistent; cannot load");
          return;
        }

      }
      //{
      //  auto temp = cbica::normPath(m_tempFolderLocation + "/temp_prev.nii.gz");
      //  cbica::WriteImage< ImageTypeFloat3D >(mSlicerManagers[0]->mITKImage, temp);
      //  auto imageInfoPrev = cbica::ImageInfo(temp);
      //  auto sizePrev = imageInfoPrev.GetImageSize();
      //  auto size = imageInfo.GetImageSize();
      //  const std::string errorMsg = " not matching. Please register the image(s). Skipping file: ";
      //  for (size_t i = 0; i < ImageTypeFloat3D::ImageDimension; i++)
      //  {
      //    if (sizePrev[i] != size[i])
      //    {
      //      updateProgress(0, "Size" + errorMsg + fname);
      //      ShowErrorMessage("Size" + errorMsg + fname);
      //      return; //
      //    }
      //  }
      //  auto spacingPrev = imageInfoPrev.GetImageSpacings();
      //  auto spacing = imageInfo.GetImageSpacings();
      //  for (size_t i = 0; i < ImageTypeFloat3D::ImageDimension; i++)
      //  {
      //    if (spacing[i] != spacingPrev[i])
      //    {
      //      updateProgress(0, "Spacing" + errorMsg + fname);
      //      ShowErrorMessage("Spacing" + errorMsg + fname);
      //      return; //
      //    }
      //  }
      //  auto originPrev = imageInfoPrev.GetImageOrigins();
      //  auto origin = imageInfo.GetImageOrigins();
      //  if (!m_advancedVisualizer)
      //  {
      //    for (size_t i = 0; i < ImageTypeFloat3D::ImageDimension; i++)
      //    {
      //      if (origin[i] != originPrev[i])
      //      {
      //        updateProgress(0, "Origin" + errorMsg + fname);
      //        ShowErrorMessage("Origin" + errorMsg + fname);
      //        return; //
      //      }
      //    }
      //  }
      //}

      if (bSkipDup)
      {
        for (int j = 0; j < (int)mSlicerManagers.size(); j++)
        {
          if (fname == mSlicerManagers[j]->GetPathFileName())
          {
            updateProgress(0, "Duplicate file skipped :" + fname);
            return;
          }
        }
      }
    }

    QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));

    if (imageInfo.GetImageDimensions() == 4)
    {
      image4DSlider->setEnabled(true);
      image4DSlider->setRange(0, imageInfo.GetImageSize()[3] - 1);
      ImageTypeFloat4D::Pointer imagePerf = cbica::ReadImage<ImageTypeFloat4D>(fname);
      imageManager->SetPerfImage(imagePerf);
      imageManager->mImageSubType = CAPTK::ImageModalityType::IMAGE_TYPE_PERFUSION;
      //return;
    }
    else
    {
      imageManager->SetOriginalOrigin(imageInfo.GetImageOrigins());
      auto currentImage = cbica::ReadImage<ImageTypeFloat3D>(fname);
      imageManager->SetOriginalDirection(currentImage->GetDirection());
      imageManager->SetOriginalOrigin(currentImage->GetOrigin());
      currentImage = ChangeImageDirectionToIdentity< ImageTypeFloat3D >(cbica::ReadImageWithOrientFix< ImageTypeFloat3D >(fname));
      imageManager->SetImage(currentImage);
      imageManager->mImageSubType = guessImageType(fname);
    }
    mInputPathName = cbica::getFilenamePath(fname).c_str();
    imageManager->SetFilename(fname);
    imageManager->SetMask(mMask);
    imageManager->setTempFolderLocation(m_tempFolderLocation);
    int rowIndex = (int)mSlicerManagers.size();

    m_imagesTable->setRowCount(rowIndex + 1);
    mSlicerManagers.push_back(imageManager);


    QFileInfo fileinfo(imageManager->GetFileName().c_str());
    QString id = fname.c_str() + QString::number(mSlicerManagers.size() - 1);
    //
    std::string strImageType = " IMAGE ";
    {
      QTableWidgetItem *item = new QTableWidgetItem(fileinfo.fileName());
      item->setData(Qt::UserRole, id.toStdString().c_str());
      item->setFlags(item->flags() & ~Qt::ItemIsEditable);

      QTablePushButton* cButton = new QTablePushButton;
      cButton->setItem(item);
      cButton->setText(QString("X"));
      connect(cButton, SIGNAL(clickedInto(QTableWidgetItem*)), this, SLOT(CloseImage(QTableWidgetItem*)));

      QLabel * label = new QLabel;
      label->setText(QString::fromStdString(strImageType));
      m_imagesTable->setCellWidget(rowIndex, TAB_IMAGES_COLUMN_CLOSE, cButton);
      m_imagesTable->setCellWidget(rowIndex, TAB_IMAGES_COLUMN_TYPE, label);
      m_imagesTable->setItem(rowIndex, TAB_IMAGES_COLUMN_NAME, item);
    }

    imagesPanel->NewImageLoaded(id, imageManager->GetBaseFileName(), rowIndex, strImageType, imageManager->mImageSubType, this);



    QTableWidgetItem *item = new QTableWidgetItem(fileinfo.fileName());
    item->setData(Qt::UserRole, id.toStdString().c_str());

    QTablePushButton* cButton = new QTablePushButton;
    cButton->setItem(item);
    cButton->setText(QString("X"));
    connect(cButton, SIGNAL(clickedInto(QTableWidgetItem*)), this, SLOT(CloseImage(QTableWidgetItem*)));

    mSlicerManagers.back()->SetId(id.toStdString());
    connect(mSlicerManagers.back(), SIGNAL(LeftButtonReleaseSignal(int)), this, SLOT(propogateSlicerPosition(int)));
    connect(mSlicerManagers.back(), SIGNAL(currentImageChanged(std::string &)), this, SLOT(CurrentImageChanged(std::string &)));
    connect(mSlicerManagers.back(), SIGNAL(currentPickedImageChanged(std::string)), this, SLOT(CurrentPickedImageChanged(std::string)));
    connect(mSlicerManagers.back(), SIGNAL(UpdatePosition(int, double, double, double, double, double, double, double)), this, SLOT(MousePositionChanged(int, double, double, double, double, double, double, double)));
    connect(mSlicerManagers.back(), SIGNAL(WindowLevelChanged()), this, SLOT(WindowLevelChanged()));
    connect(mSlicerManagers.back(), SIGNAL(UpdateSlice(int, int)), this, SLOT(UpdateSlice(int, int)));
    connect(mSlicerManagers.back(), SIGNAL(UpdateSliceRange(int, int, int)), this, SLOT(UpdateSliceRange(int, int, int)));
    connect(mSlicerManagers.back(), SIGNAL(UpdateLinkManager(std::string, int, double, double, double)), this, SLOT(UpdateLinkManager(std::string, int, double, double, double)));
    connect(mSlicerManagers.back(), SIGNAL(ChangeImageWithOrder(SlicerManager*, int)), this, SLOT(ChangeImageWithOrder(SlicerManager*, int)));
    connect(mSlicerManagers.back(), SIGNAL(UpdateBorderWidgetInMain(double, double, double, double)), this, SLOT(UpdateBorderWidget(double, double, double, double)));
    connect(mSlicerManagers.back(), SIGNAL(UpdateBorderWidgetInMain(double, double)), this, SLOT(UpdateBorderWidget(double, double)));
    connect(mSlicerManagers.back(), SIGNAL(UpdateActionInMain(const QVariantList&)), this, SLOT(UpdateActionQ(const QVariantList&)));

    //connect(mSlicerManagers.back(), SIGNAL(SeedPointsAdded()), tumorPanel, SLOT(sAddPoint()));
    //connect(mSlicerManagers.back(), SIGNAL(SeedPointsAdded(int, bool)), tumorPanel, SLOT(sAddPoint(int, bool)));
    //connect(mSlicerManagers.back(), SIGNAL(TissuePointsAdded(int)), tumorPanel, SLOT(tAddPoint(int)));
    //connect(m_tabWidget, SIGNAL(currentChanged(int)), tumorPanel, SLOT(tabSelected()));
    InitSlicers();

    if (bFirstLoad)
    {
      InitMask(mSlicerManagers.back()->mImage);
    }
    for (int j = 0; j < (int)mSlicerManagers.back()->mSlicers.size(); j++)
    {
      mSlicerManagers.back()->mSlicers[j]->SetMask(mSlicerManagers.back()->GetMask());
    }
    if (fname.find("scan_label_map") != std::string::npos)
    {
      mSlicerManagers.back()->SetPreset(PRESET_LABEL);
    }
    if (fname.find("gt") != std::string::npos)
    {
      mSlicerManagers.back()->SetPreset(PRESET_LABEL2);
    }
    if (fname.find("roiDE") != std::string::npos)
    {
      mSlicerManagers.back()->SetPreset(PRESET_PROB);
    }

    if (mSlicerManagers.size() > 0)
    {
      if (mSlicerManagers.back()->mMask->GetDimensions()[2] != 1)
      {
        CoronalViewWidget->show();
        SaggitalViewWidget->show();
      }
      AxialViewWidget->show();
      infoPanel->show();

      windowLabel->setEnabled(true);
      windowSpinBox->setEnabled(true);
      levelLabel->setEnabled(true);
      levelSpinBox->setEnabled(true);
      presetLabel->setEnabled(true);
      presetComboBox->setEnabled(true);

      if (bFirstLoad)
      {
        for (int i = 0; i < 3; i++)
        {
          mSlicerManagers.back()->GetSlicer(i)->SetInitPosition();
        }
        DisplayChanged(item);
      }
      else
      {
        QTableWidgetItem* item = NULL;
        for (int i = 0; i < (int)mSlicerManagers.size(); i++)
        {
          item = GetItemFromSlicerManager(mSlicerManagers[i]);
          if (!item->isSelected())
          {
            item->setSelected(true);
          }
        }
        DisplayChanged(item);
      }

      if (mSlicerManagers.size() > 1)
      {
        for (int i = 0; i < (int)mSlicerManagers.size(); i++)
        {
          for (int j = i + 1; j < (int)mSlicerManagers.size(); j++)
          {
            AddLink(/*QString::fromStdString*/(mSlicerManagers[i]->GetId().c_str()), /*QString::fromStdString*/(mSlicerManagers[j]->GetId().c_str()));
          }
        }
      }
      QTableWidgetItem* item = GetItemFromSlicerManager(mSlicerManagers.back());
      item->setSelected(true);
      InitDisplay();
    }
    propogateSlicerPosition();
    updateProgress(0);
    QApplication::restoreOverrideCursor();
  }
//  else
//  {
//    auto path = cbica::getFilenamePath(fileName);
//    dicomfilename = fileName;
//    auto filesInDir = cbica::filesInDirectory(path, false);
//
//#ifndef _WIN32
//    for (auto it = filesInDir.begin(); it != filesInDir.end();)
//    {
//      if ((*it == ".") || (*it == ".."))
//      {
//        filesInDir.erase(it);
//      }
//      else
//      {
//        ++it;
//      }
//    }
//#endif
//
//    // remove any files that aren't DICOM (thumbs.db and stuff like that)
//    for (size_t i = 0; i < filesInDir.size(); i++)
//    {
//      if (cbica::getFilenameExtension(path + "/" + filesInDir[i]) != ".dcm")
//      {
//        filesInDir.erase(filesInDir.begin() + i);
//      }
//    }
//
//    if (filesInDir.size() == 1) // single DICOM slice
//    {
//      ConversionFrom2Dto3D(fileName, true);
//    }
//    else // for 3D images, call dcm2nii
//    {
//      CallDCM2NIfTIConversion(fileName, true);
//    }
//
//    return;
//  }
}

void fMainWindow::CurrentImageChanged(std::string &id)
{
  if (id == mCurrentSelectedImageId)
  {
    return;
  }
  int selected = 0;
  for (int i = 0; i < m_imagesTable->rowCount(); i++)
  {
    if (m_imagesTable->item(i, TAB_IMAGES_COLUMN_NAME)->data(Qt::UserRole).toString().toStdString() == id)
    {
      selected = i;
    }
    else
    {
      m_imagesTable->item(i, TAB_IMAGES_COLUMN_NAME)->setSelected(false);
    }
  }
  m_imagesTable->item(selected, TAB_IMAGES_COLUMN_NAME)->setSelected(true);
  mCurrentSelectedImageId = id;
  emit SelectedImageHasChanged(mSlicerManagers[selected]);
  m_imagesTable->resizeColumnsToContents();
  m_imagesTable->resizeRowsToContents();
}

void fMainWindow::propogateSlicerPosition(int slicerId, int imageId)
{
  //TBD this not a proper fix a lot work needs to be done to make slicer, slicerManager slicerManagerCommand to be made clean OOP
  if (imageId < 0)
  {
    auto items = m_imagesTable->selectedItems();
    if (!items.empty())
    {
      imageId = GetSlicerIndexFromItem(items[0]);
    }

  }
  if (imageId < 0 || imageId >= (int)mSlicerManagers.size())
  {
    return;
  }
  const int MAX_SLICES = 3;
  if (slicerId < 0 || slicerId >= MAX_SLICES)
  {
    return;
  }
  double* pos = mSlicerManagers[imageId]->GetSlicer(slicerId)->GetCurrentPosition();
  for (size_t r = 0; r < mSlicerManagers.size(); r++)
  {
    for (size_t i = 0; i < MAX_SLICES; i++)
    {
      mSlicerManagers[r]->GetSlicer(i)->SetCurrentPosition(pos[0], pos[1], pos[2]);
    }
  }
}
void fMainWindow::CurrentPickedImageChanged(std::string id)
{
  if (id == mCurrentPickedImageId) {
    return;
  }
  int selected = 0;
  for (int i = 0; i < m_imagesTable->rowCount(); i++) {
    if (m_imagesTable->item(i, TAB_IMAGES_COLUMN_NAME)->data(Qt::UserRole).toString().toStdString() == id) {
      selected = i;
    }
    else {
      m_imagesTable->item(i, TAB_IMAGES_COLUMN_NAME)->setSelected(false);
    }
  }
  m_imagesTable->item(selected, TAB_IMAGES_COLUMN_NAME)->setSelected(true);
  mCurrentPickedImageId = id;
  mCurrentPickedImageIndex = selected;
}

void fMainWindow::ImageInfoChanged()
{
  auto items = m_imagesTable->selectedItems();
  if (items.empty()) {
    return;
  }
  int index = GetSlicerIndexFromItem(items[0]);
  if (index < 0 || index >= (int)mSlicerManagers.size()) {
    return;
  }

  vtkSmartPointer<vtkMatrix4x4> transformation;
  QString image = m_imagesTable->selectedItems()[0]->data(Qt::DisplayRole).toString();

  vtkSmartPointer<vtkImageData> imageSelected;
  vtkSmartPointer<vtkTransform> transformSelected;
  imageSelected = mSlicerManagers[index]->GetSlicer(0)->GetImage();
  transformSelected = mSlicerManagers[index]->GetSlicer(0)->GetTransform();

  std::string vtktype = mSlicerManagers[index]->GetImage()->GetScalarTypeAsString();
  QString pixelType = vtktype.c_str();

  infoPanel->setFileName(image);
  infoPanel->setSizePixel(imageSelected->GetDimensions()[0], imageSelected->GetDimensions()[1], imageSelected->GetDimensions()[2]);
  infoPanel->setOrigin(mSlicerManagers[index]->mOrigin[0], mSlicerManagers[index]->mOrigin[1], mSlicerManagers[index]->mOrigin[2]);
  infoPanel->setSpacing(imageSelected->GetSpacing()[0], imageSelected->GetSpacing()[1], imageSelected->GetSpacing()[2]);
  transformation = transformSelected->GetMatrix();
  //tumorPanel->SetCurrentSPoints(mSlicerManagers[index]->mSeedPoints);
  //tumorPanel->SetCurrentTPoints(mSlicerManagers[index]->mTissuePoints);
  //tumorPanel->SetCurrentPath(mInputPathName.toStdString());

  for (int i = 0; i < 3; i++)
  {
    mSlicerManagers[index]->UpdateInfoOnCursorPosition(i);
  }

  WindowLevelChanged();

  // reset SliceManager order
  for (int j = 0; j < (int)mSlicerManagers.size(); j++) {
    mSlicerManagers[j]->SetOrder(j);
  }
  for (int i = 0; i < m_imagesTable->rowCount(); i++) {
    QString id_table = m_imagesTable->item(i, TAB_IMAGES_COLUMN_NAME)->data(Qt::UserRole).toString();
    for (int j = 0; j < (int)mSlicerManagers.size(); j++) {
      QString id_sm = mSlicerManagers[j]->GetId().c_str();
      if (id_table == id_sm) {
        mSlicerManagers[j]->SetOrder(i);
        break;
      }
    }
  }
}

void fMainWindow::DisplayChanged()
{
  auto items = m_imagesTable->selectedItems();
  if (items.empty())
  {
    return;
  }
  DisplayChanged(items[0]);
}
void fMainWindow::DisplayChanged(QTableWidgetItem *clickedItem)
{
  int slicerManagerIndex = GetSlicerIndexFromItem(clickedItem);
  if (slicerManagerIndex < 0 || slicerManagerIndex >= (int)mSlicerManagers.size())
  {
    return;
  }

  QTableWidgetItem* clickedParentItem = m_imagesTable->item(slicerManagerIndex, TAB_IMAGES_COLUMN_NAME);
  for (int i = 0; i < (int)mSlicerManagers.size(); i++)
  {
    QTableWidgetItem* currentParentItem = m_imagesTable->item(i, TAB_IMAGES_COLUMN_NAME);
    if (currentParentItem != clickedParentItem)
    {
      for (int j = 0; j < 3; j++)
      {
        mSlicerManagers[i]->UpdateSlicer(j, false);
      }
    }
    else
    {
      int VisibleInWindow = 0;
      for (int j = 0; j < 3; j++)
      {
        mSlicerManagers[i]->UpdateSlicer(j, true);
        mSlicerManagers[i]->UpdateInfoOnCursorPosition(j);
        DisplaySliders(i, j);
        //
        if (mSlicerManagers[i]->GetSlicer(j)->GetActive())
        {
          VisibleInWindow = j;
        }
      }

      mSlicerManagers[i]->Picked();
      mSlicerManagers[i]->UpdateViews(VisibleInWindow);
      mSlicerManagers[i]->UpdateLinked(VisibleInWindow);
      mSlicerManagers[i]->UpdateInfoOnCursorPosition(VisibleInWindow);
    }

    for (int j = 0; j < 3; j++)
    {
      mSlicerManagers[i]->mSlicers[j]->RemoveOverlay();
    }
  }
  UpdateRenderWindows();

  ImageInfoChanged();
}



int fMainWindow::GetSlicerIndexFromItem(QTableWidgetItem* item)
{
  if (item != NULL) {
    QString id = item->data(Qt::UserRole).toString();
    std::string id_string = id.toStdString();
    for (int i = 0; i < m_imagesTable->rowCount(); i++)
    {
      QString id_table_string = m_imagesTable->item(i, TAB_IMAGES_COLUMN_NAME)->data(Qt::UserRole).toString();
      std::string table_string = id_table_string.toStdString();
      if (m_imagesTable->item(i, TAB_IMAGES_COLUMN_NAME)->data(Qt::UserRole).toString() == id)
      {
        return i;
      }
    }
  }
  return -1;
}

QTableWidgetItem* fMainWindow::GetItemFromSlicerManager(SlicerManager* sm)
{
  QString id = sm->GetId().c_str();
  for (int i = 0; i < m_imagesTable->rowCount(); i++)
  {
    if (m_imagesTable->item(i, TAB_IMAGES_COLUMN_NAME)->data(Qt::UserRole).toString() == id)
    {
      return m_imagesTable->item(i, TAB_IMAGES_COLUMN_NAME);
    }
  }
  return NULL;
}

void fMainWindow::InitSlicers()
{
  if (mSlicerManagers.size()) {
    mSlicerManagers.back()->GenerateDefaultLookupTable();

    mSlicerManagers.back()->SetSlicerWindow(0, AxialViewWidget->GetRenderWindow());
    mSlicerManagers.back()->SetSlicerWindow(1, CoronalViewWidget->GetRenderWindow());
    mSlicerManagers.back()->SetSlicerWindow(2, SaggitalViewWidget->GetRenderWindow());
  }
}

void fMainWindow::InitDisplay()
{
  if (!mSlicerManagers.empty())
  {
    for (int j = 0; j < 3; j++)
    {
      InteractorStyleNavigator* style = InteractorStyleNavigator::New();
      style->SetAutoAdjustCameraClippingRange(1);
      for (int i = 0; i < m_imagesTable->rowCount(); i++)
      {
        mSlicerManagers[i]->SetInteractorStyleNavigator(j, style);
        //
        mSlicerManagers[i]->updateToRefCam(mSlicerManagers[i]->GetSlicer(0));
        mSlicerManagers[i]->GetSlicer(j)->SetInitPosition();
      }
      style->Delete();
    }
  }
}

void fMainWindow::DisplaySliders(int slicer, int window)
{
  int range[2];
  mSlicerManagers[slicer]->GetSlicer(window)->GetSliceRange(range);
  int position = mSlicerManagers[slicer]->GetSlicer(window)->GetSlice();

  bool showVertical = false;
  if (GetNumberOfDimensions(mSlicerManagers[slicer]->GetSlicer(window)->GetImage()) >= 3) {
    showVertical = true;
  }
  if (showVertical) {
    verticalSliders[window]->show();
  }
  else {
    verticalSliders[window]->hide();
  }
  verticalSliders[window]->setRange(range[0], range[1]);
  verticalSliders[window]->setValue(position);
}

void fMainWindow::CloseImage(QTableWidgetItem* item)
{
  int index = GetSlicerIndexFromItem(item);
  if (index < 0 || index >= (int)mSlicerManagers.size())
  {
    return;
  }
  if (mSlicerManagers.size() > 1)
  {
    for (int k = 0; k < (int)mSlicerManagers.size() - 1; k++)
    {
      if (k != index)
      {
        RemoveLink(/*QString::fromStdString*/(mSlicerManagers[k]->GetId()), /*QString::fromStdString*/(mSlicerManagers[index]->GetId()));
      }
    }
  }
  std::vector<SlicerManager*>::iterator Manageriter = mSlicerManagers.begin();
  for (int i = 0; i < index; i++)
  {
    Manageriter++;
  }

  mSlicerManagers[index]->RemoveActors();
  mSlicerManagers[index]->mTissuePoints->Clear();
  mSlicerManagers[index]->mSeedPoints->Clear();
  delete mSlicerManagers[index];
  mSlicerManagers.erase(Manageriter);
  m_imagesTable->removeRow(index);

  if (mSlicerManagers.size() >= 1)
  {
    QTableWidgetItem* item_tmp = GetItemFromSlicerManager(mSlicerManagers.back());
    item_tmp->setSelected(true);
    DisplayChanged(item_tmp);
  }
  else
  {
    this->ResetNumberOfPoints();
    AxialViewWidget->hide();
    CoronalViewWidget->hide();
    SaggitalViewWidget->hide();

    for (int i = 0; i < 3; i++)
    {
      verticalSliders[i]->hide();
    }
    //
    windowLabel->setEnabled(false);
    windowSpinBox->setEnabled(false);
    levelLabel->setEnabled(false);
    levelSpinBox->setEnabled(false);
    thresholdLabel->setEnabled(false);
    thresholdSpinBox->setEnabled(false);
    presetLabel->setEnabled(false);
    presetComboBox->setEnabled(false);
    m_imgGeodesicOut = NULL;
    mLandmarks->Clear();
    mSeedPoints->Clear();
    mTissuePoints->Clear();
  }

  InitDisplay();
}

void fMainWindow::MousePositionChanged(int visibility, double x, double y, double z, double X, double Y, double Z, double value)
{
  infoPanel->setCurrentInfo(visibility, x, y, z, X, Y, Z, value);
  //tumorPanel->HighlightCurrentSelctedPoints(x, y, z, X, Y, Z, value);
}

void fMainWindow::WindowLevelChanged()
{
  auto items = m_imagesTable->selectedItems();
  if (items.empty()) {
    return;
  }
  int index = GetSlicerIndexFromItem(items[0]);
  if (index < 0 || index >= (int)mSlicerManagers.size())
  {
    return;
  }
  windowSpinBox->setValue(mSlicerManagers[index]->GetColorWindow());
  levelSpinBox->setValue(mSlicerManagers[index]->GetColorLevel());
  thresholdSpinBox->setValue(mSlicerManagers[index]->GetThresholdIndex());
  presetComboBox->setCurrentIndex(mSlicerManagers[index]->GetPreset());
  if (presetComboBox->currentIndex() == PRESET_THRESHOLD || presetComboBox->currentIndex() == PRESET_GEODESIC)
  {
    thresholdLabel->setEnabled(true);
    thresholdSpinBox->setEnabled(true);
  }
  else
  {
    thresholdLabel->setEnabled(false);
    thresholdSpinBox->setEnabled(false);
  }
}

void fMainWindow::WindowLevelEdited()
{
  presetComboBox->setCurrentIndex(PRESET_USER);
  UpdateWindowLevel();
}

void fMainWindow::SetWindowLevel(double w, double l)
{
  windowSpinBox->setValue(w);
  levelSpinBox->setValue(l);
  presetComboBox->setCurrentIndex(PRESET_USER);
  UpdateWindowLevel();
}

void fMainWindow::UpdateWindowLevel()
{
  if (!m_ComparisonMode)//! if comparison mode OFF
  {
    auto items = m_imagesTable->selectedItems();
    if (items.empty()) {
      return;
    }
    int index = GetSlicerIndexFromItem(items[0]);
    if (index >= 0 && index < (int)mSlicerManagers.size())
    {
      mSlicerManagers[index]->SetColorWindow(windowSpinBox->value());
      mSlicerManagers[index]->SetColorLevel(levelSpinBox->value());
      mSlicerManagers[index]->SetPreset(presetComboBox->currentIndex());
      mSlicerManagers[index]->Render();
      //
      if (presetComboBox->currentIndex() == PRESET_THRESHOLD || presetComboBox->currentIndex() == PRESET_GEODESIC) {
        thresholdLabel->setEnabled(true);
        thresholdSpinBox->setEnabled(true);
      }
      else
      {
        thresholdLabel->setEnabled(false);
        thresholdSpinBox->setEnabled(false);
      }
      //
      WindowLevelChanged();
    }
  }
  else//! if comparison mode ON
  {
    std::vector<vtkSmartPointer<Slicer>> comparisonViewers = this->GetComparisonViewers();
    for (int i = 0; i < comparisonViewers.size(); i++)
    {
      comparisonViewers[i]->SetColorWindow(windowSpinBox->value());
      comparisonViewers[i]->SetColorLevel(levelSpinBox->value());
      comparisonViewers[i]->Render();
    }
  }
}

void fMainWindow::thresholdSpinBoxChanged()
{
  if (presetComboBox->currentIndex() == PRESET_THRESHOLD)
  {
    auto items = m_imagesTable->selectedItems();
    if (items.empty())
    {
      return;
    }
    int index = GetSlicerIndexFromItem(items[0]);
    if (index >= 0 && index < (int)mSlicerManagers.size())
    {
      mSlicerManagers[index]->SetThresholdIndex(thresholdSpinBox->value());
      mSlicerManagers[index]->SetPreset(mSlicerManagers[index]->GetPreset());
      mSlicerManagers[index]->Render();
      WindowLevelChanged();
    }
  }
  else if (presetComboBox->currentIndex() == PRESET_GEODESIC)
  {
    ApplicationGeodesicTreshold();
  }

}

void fMainWindow::UpdateLinkManager(std::string id, int slicer, double x, double y, double z)
{
  for (int i = 0; i < (int)mSlicerManagers.size(); i++)
  {
    if (mSlicerManagers[i]->GetId() == id)
    {
      mSlicerManagers[i]->GetSlicer(slicer)->SetCurrentPosition(x, y, z);
      mSlicerManagers[i]->UpdateViews(slicer);
      break;
    }
  }
}
void fMainWindow::UpdateLinkedNavigation(Slicer* refSlicer)
{
  for (int i = 0; i < (int)mSlicerManagers.size(); i++)
  {
    mSlicerManagers[i]->updateToRefCam(refSlicer);
  }
}

void fMainWindow::CloseImage()
{
  auto items = m_imagesTable->selectedItems();
  if (items.empty()) {
    return;
  }
  CloseImage(items[0]);
}
void fMainWindow::CloseAllImages()
{
  while (m_imagesTable->rowCount() > 0)
  {
    m_imagesTable->selectRow(0);//TBD speedup this
    CloseImage();
  }
  infoPanel->setFileName("");
  infoPanel->setSizePixel(0, 0, 0);
  infoPanel->setOrigin(0, 0, 0);
  infoPanel->setSpacing(0, 0, 0);
  infoPanel->setCurrentInfo(0, 0, 0, 0, 0, 0, 0, 0);
}
void fMainWindow::ResetTransformationToIdentity()
{
  auto items = m_imagesTable->selectedItems();
  if (items.empty()) {
    return;
  }
  int index = GetSlicerIndexFromItem(items[0]);
  if (index >= 0 && index < (int)mSlicerManagers.size()) {
    mSlicerManagers[index]->ResetTransformationToIdentity();
    ImageInfoChanged();
  }
}

void fMainWindow::AddLink(const std::string &image1, const std::string &image2)
{

  for (int i = 0; i < (int)mSlicerManagers.size(); i++)
  {
    if (image1/*.toStdString()*/ == mSlicerManagers[i]->GetId())
    {
      mSlicerManagers[i]->AddLink(image2/*.toStdString()*/);
      //sm1 = i;
    }
    if (image2/*.toStdString()*/ == mSlicerManagers[i]->GetId())
    {
      mSlicerManagers[i]->AddLink(image1/*.toStdString()*/);
      //sm2 = i;
    }
  }
}

void fMainWindow::RemoveLink(const std::string &image1, const std::string &image2)
{
  for (int i = 0; i < (int)mSlicerManagers.size(); i++)
  {
    if (image1/*.toStdString()*/ == mSlicerManagers[i]->GetId())
    {
      mSlicerManagers[i]->RemoveLink(image2/*.toStdString()*/);
    }
    if (image2/*.toStdString()*/ == mSlicerManagers[i]->GetId()) {
      mSlicerManagers[i]->RemoveLink(image1/*.toStdString()*/);
    }
  }
}

void fMainWindow::ChangeImageWithOrder(SlicerManager *sm, int order)
{
  if (mSlicerManagers.size() <= 1) {
    return;
  }
  if (order >= (int)mSlicerManagers.size()) {
    return;
  }

  QTableWidgetItem* item;
  item = GetItemFromSlicerManager(mSlicerManagers[order]);
  item->setSelected(true);
  DisplayChanged(item);
}

void fMainWindow::SetImageInfoIntensityValue(double value)
{
  this->infoPanel->setIntensityValue(value);
}

void fMainWindow::SetImageInfoZSlicePosition(int zslice)
{
  this->infoPanel->setZSlicePosition(zslice);
}

void fMainWindow::OnSliderMovedInComparisonMode(int value)
{
  if (AxialViewSlider->value() != value || 
    CoronalViewSlider->value() != value ||
    SaggitalViewSlider->value() != value)
  {
    AxialViewSlider->setValue(value);
    CoronalViewSlider->setValue(value);
    SaggitalViewSlider->setValue(value);

  }
  if (m_ComparisonViewerLeft->GetSlice() != value)
  {
    m_ComparisonViewerLeft->SetSlice(value);
    m_ComparisonViewerCenter->SetSlice(value);
    m_ComparisonViewerRight->SetSlice(value);

    m_ComparisonViewerLeft->Render();
    m_ComparisonViewerCenter->Render();
    m_ComparisonViewerRight->Render();
  }

}

void fMainWindow::AxialViewSliderChanged()
{
  static int value = -1;

  value = AxialViewSlider->value();
  auto items = m_imagesTable->selectedItems();
  if (items.empty()) {
    return;
  }
  int index = GetSlicerIndexFromItem(items[0]);
  if (index >= 0 && index < (int)mSlicerManagers.size()) {
    if (mSlicerManagers[index]->GetSlicer(0)->GetSlice() != value) {
      mSlicerManagers[index]->GetSlicer(0)->SetSlice(value);
      mSlicerManagers[index]->VerticalSliderHasChanged(0, value);
      mSlicerManagers[index]->UpdateSlice(0);
    }
  }
}

void fMainWindow::SaggitalViewSliderChanged()
{
  static int value = -1;
  if (value == SaggitalViewSlider->value())
  {
    return;
  }
  else {
    value = SaggitalViewSlider->value();
  }
  auto items = m_imagesTable->selectedItems();
  if (items.empty()) {
    return;
  }
  int index = GetSlicerIndexFromItem(items[0]);
  if (index >= 0 && index < (int)mSlicerManagers.size()) {
    if (mSlicerManagers[index]->GetSlicer(2)->GetSlice() != value) {
      mSlicerManagers[index]->GetSlicer(2)->SetSlice(value);
      mSlicerManagers[index]->VerticalSliderHasChanged(2, value);
      mSlicerManagers[index]->UpdateSlice(2);
    }
  }
}

void fMainWindow::CoronalViewSliderChanged()
{
  static int value = -1;
  if (value == CoronalViewSlider->value())
  {
    return;
  }
  else
  {
    value = CoronalViewSlider->value();
  }
  auto items = m_imagesTable->selectedItems();
  if (items.empty()) {
    return;
  }
  int index = GetSlicerIndexFromItem(items[0]);
  if (index >= 0 && index < (int)mSlicerManagers.size()) {
    if (mSlicerManagers[index]->GetSlicer(1)->GetSlice() != value) {
      mSlicerManagers[index]->GetSlicer(1)->SetSlice(value);
      mSlicerManagers[index]->VerticalSliderHasChanged(1, value);
      mSlicerManagers[index]->UpdateSlice(1);
    }
  }
}

void fMainWindow::UpdateSlice(int slicer, int slice)
{
  if (slicer == 0)
  {
    AxialViewSlider->setValue(slice);
  }
  else if (slicer == 1)
  {
    CoronalViewSlider->setValue(slice);
  }
  else if (slicer == 2)
  {
    SaggitalViewSlider->setValue(slice);
  }
  propogateSlicerPosition();
}

void fMainWindow::UpdateSliceRange(int slicer, int min, int max)
{
  int position = int((min + max) / 2);
  if (slicer == 0) {
    AxialViewSlider->setValue(position);
    AxialViewSlider->setRange(min, max);
  }
  else if (slicer == 1) {
    CoronalViewSlider->setValue(position);
    CoronalViewSlider->setRange(min, max);
  }
  else if (slicer == 2) {
    SaggitalViewSlider->setValue(position);
    SaggitalViewSlider->setRange(min, max);
  }
}

void fMainWindow::UpdateRenderWindows()
{
  if (m_imagesTable->rowCount() <= 0)
  {
    return;
  }
  auto items = m_imagesTable->selectedItems();
  if (items.empty()) {
    return;
  }
  QTableWidgetItem* item = items[0];
  if (item == NULL) {
    return;
  }
  int index = GetSlicerIndexFromItem(item);
  if (index >= 0 && index < (int)mSlicerManagers.size())
  {
    mSlicerManagers[index]->Render();
  }
  //*/
}

void fMainWindow::SetActiveLandmarksType(int type, int row, int col)
{
  mCurrentLandmarkTissueType = row;
  if (type == LANDMARK_TYPE::DEFAULT)
  {
    emit LandmarksFocused(true);
    emit SeedPointsFocused(false);
    emit TissuePointsFocused(false);
  }
  else if (type == LANDMARK_TYPE::TUMOR_POINTS)
  {
    emit LandmarksFocused(false);
    emit SeedPointsFocused(true);
    emit TissuePointsFocused(false);
    //tumorPanel->mTumorPointsSelected = true;
  }
  else if (type == LANDMARK_TYPE::TISSUE_POINTS)
  {
    emit LandmarksFocused(false);
    emit SeedPointsFocused(false);
    emit TissuePointsFocused(true);
  }
  else
  {
    emit LandmarksFocused(false);
    emit SeedPointsFocused(false);
    emit TissuePointsFocused(false);
  }

  for (int i = 0; i < (int)mSlicerManagers.size(); i++)
    mSlicerManagers[i]->SetCurrentLandmarksType(type, row, col);

  UpdateRenderWindows();
}
void fMainWindow::panelChanged(int current)
{
  if (drawingPanel) //Reset shape mode on every panel switch
  {
    m_drawShapeMode = SHAPE_MODE_NONE;
    drawingPanel->shapesNoneButtonFunctionality();
  }

  else
  {
    SetActiveLandmarksType(LANDMARK_TYPE::NONE, 0, 0);
  }
  //else if (current == TAB_TUMOR)
  //{
  //  SetActiveLandmarksType(LANDMARK_TYPE::NONE, 0, 0);
  //  //tumorPanel->SetCurrentSelectedTissueType();
  //}
}

void fMainWindow::MoveSlicerCursor(double x, double y, double z, int mode)
{
  if (mCurrentPickedImageIndex < 0 || mCurrentPickedImageIndex >= (int)mSlicerManagers.size()) {
    return;
  }
  //
  if (mode == 0)
  {
    // x, y, z are LPS
    mSlicerManagers[mCurrentPickedImageIndex]->GetSlicer(0)->SetCurrentPosition(x, y, z);
    //
    mSlicerManagers[mCurrentPickedImageIndex]->Picked();
    mSlicerManagers[mCurrentPickedImageIndex]->UpdateViews(0);
    mSlicerManagers[mCurrentPickedImageIndex]->UpdateLinked(0);
    mSlicerManagers[mCurrentPickedImageIndex]->UpdateInfoOnCursorPosition(0);
  }
  else if (mode == 1)
  {
    // x, y, z are pixel
    x = x * mSlicerManagers[mCurrentPickedImageIndex]->GetSlicer(0)->GetInput()->GetSpacing()[0] + mSlicerManagers[mCurrentPickedImageIndex]->GetSlicer(0)->GetInput()->GetOrigin()[0];
    y = y * mSlicerManagers[mCurrentPickedImageIndex]->GetSlicer(0)->GetInput()->GetSpacing()[1] + mSlicerManagers[mCurrentPickedImageIndex]->GetSlicer(0)->GetInput()->GetOrigin()[1];
    z = z * mSlicerManagers[mCurrentPickedImageIndex]->GetSlicer(0)->GetInput()->GetSpacing()[2] + mSlicerManagers[mCurrentPickedImageIndex]->GetSlicer(0)->GetInput()->GetOrigin()[2];
    //
    mSlicerManagers[mCurrentPickedImageIndex]->GetSlicer(0)->SetCurrentPosition(x, y, z);
    //
    mSlicerManagers[mCurrentPickedImageIndex]->Picked();
    mSlicerManagers[mCurrentPickedImageIndex]->UpdateViews(0);
    mSlicerManagers[mCurrentPickedImageIndex]->UpdateLinked(0);
    mSlicerManagers[mCurrentPickedImageIndex]->UpdateInfoOnCursorPosition(0);
  }
  propogateSlicerPosition();
}

void fMainWindow::toolTabDockChanged(bool bUnDocked)
{
  if (bUnDocked)
  {
    m_tabWidget->setMaximumHeight(m_tabWidget->minimumHeight() * 10);
    m_toolTabdock->show();
  }
  else
  {
    m_tabWidget->setMaximumHeight(m_tabWidget->minimumHeight());
  }
}

VectorVectorDouble fMainWindow::FormulateDrawingPointsForEdemaSegmentation()
{
  VectorVectorDouble Indices;
  if (mSlicerManagers.size() <= 0)
    return Indices;

  ImageTypeShort3D::Pointer img = convertVtkToItk<short, 3>(mSlicerManagers[0]->mMask);
  typedef itk::ImageRegionIteratorWithIndex <ImageTypeShort3D> IteratorType;
  IteratorType maskIt(img, img->GetLargestPossibleRegion());
  maskIt.GoToBegin();
  while (!maskIt.IsAtEnd())
  {
    if (maskIt.Get() == 1 || maskIt.Get() == 2)
    {
      VectorDouble localIndex;
      localIndex.push_back(maskIt.GetIndex()[0]);
      localIndex.push_back(maskIt.GetIndex()[1]);
      localIndex.push_back(maskIt.GetIndex()[2]);

      Indices.push_back(localIndex);
    }
    ++maskIt;
  }
  return Indices;
}

VectorVectorDouble fMainWindow::GetMaskLabelIndices(const int label)
{
  VectorVectorDouble Indices;
  if (mSlicerManagers.size() <= 0)
    return Indices;

  ImageTypeShort3D::Pointer img = convertVtkToItk<short, 3>(mSlicerManagers[0]->mMask);
  typedef itk::ImageRegionIteratorWithIndex <ImageTypeShort3D> IteratorType;
  IteratorType maskIt(img, img->GetLargestPossibleRegion());
  maskIt.GoToBegin();
  while (!maskIt.IsAtEnd())
  {
    if (maskIt.Get() == label)
    {
      VectorDouble localIndex;
      localIndex.push_back(maskIt.GetIndex()[0]);
      localIndex.push_back(maskIt.GetIndex()[1]);
      localIndex.push_back(maskIt.GetIndex()[2]);

      Indices.push_back(localIndex);
    }
    ++maskIt;
  }
  return Indices;
}
VectorVectorDouble fMainWindow::FormulateDrawingPointsForTumorSegmentation()
{
  VectorVectorDouble Indices;
  if (mSlicerManagers.size() <= 0)
    return Indices;

  ImageTypeShort3D::Pointer img = convertVtkToItk<short, 3>(mSlicerManagers[0]->mMask);
  typedef itk::ImageRegionIteratorWithIndex <ImageTypeShort3D> IteratorType;
  IteratorType maskIt(img, img->GetLargestPossibleRegion());
  maskIt.GoToBegin();
  while (!maskIt.IsAtEnd())
  {
    if (maskIt.Get() == 1)
    {
      VectorDouble localIndex;
      localIndex.push_back(maskIt.GetIndex()[0]);
      localIndex.push_back(maskIt.GetIndex()[1]);
      localIndex.push_back(maskIt.GetIndex()[2]);

      Indices.push_back(localIndex);
    }
    ++maskIt;
  }
  return Indices;
}


void fMainWindow::SaveDicomImage()
{
  auto items = m_imagesTable->selectedItems();
  if (items.empty())
    return;
  int index = GetSlicerIndexFromItem(items[0]);

  typedef short PixelType;
  const unsigned int Dimensions = 3;

  QString directory = getExistingDirectory(this, mInputPathName);
  if (directory.isNull())
    return;

  std::string directory_wrap = directory.toStdString();
  if (directory_wrap[directory_wrap.length() - 1] != '/')
  {
    directory_wrap += "/";
  }

  typedef itk::Image<PixelType, Dimensions> ImageType;
  typedef itk::CastImageFilter< itk::Image< float, Dimensions >, ImageType > CastFilterType;
  CastFilterType::Pointer castFilter = CastFilterType::New();
  castFilter->SetInput(mSlicerManagers[index]->mITKImage);
  castFilter->Update();

  try
  {
    updateProgress(0, "Image saved! (" + directory_wrap + ")");
  }
  catch (itk::ExceptionObject & excp)
  {
    ShowErrorMessage("Couldn't write mask as DICOM: " + std::string(excp.GetDescription()));
    return;
  }

}

void fMainWindow::SaveDicomDrawing()
{
  auto items = m_imagesTable->selectedItems();
  if (items.empty())
    return;
  int index = GetSlicerIndexFromItem(items[0]);

  typedef short MaskPixelType;
  const unsigned int Dimensions = 3;

  typedef itk::Image<MaskPixelType, Dimensions> ImageTypeMask;
  ImageTypeMask::Pointer imageToWrite = convertVtkToItk<MaskPixelType, Dimensions>(mSlicerManagers[index]->mMask);

  typedef itk::MinimumMaximumImageCalculator< ImageTypeMask > CalculatorType;
  CalculatorType::Pointer calculator = CalculatorType::New();
  calculator->SetImage(imageToWrite);
  calculator->Compute();

  if (calculator->GetMaximum() == 0) // this means that the mask image contains no values at all
  {
    ShowErrorMessage("There should be at least one region to save.");
    return;
  }

  QString directory = getExistingDirectory(this, mInputPathName);
  if (directory.isNull())
    return;

  std::string directory_wrap = directory.toStdString();
  if (directory_wrap[directory_wrap.length() - 1] != '/')
  {
    directory_wrap += "/";
  }

  if (cbica::getFilenamePath(mSlicerManagers[index]->mSeriesReader->GetFileNames()[0]) == directory_wrap)
  {
    ShowErrorMessage("Cannot save to source directory. Please select another.");
    return;
  }

  try
  {
    updateProgress(0, "Success!");
  }
  catch (itk::ExceptionObject & excp)
  {
    ShowErrorMessage("Couldn't write mask as DICOM: " + std::string(excp.GetDescription()));
    return;
  }


}

void fMainWindow::SaveDrawing()
{
  auto items = m_imagesTable->selectedItems();
  if (items.empty())
    return;
  int index = GetSlicerIndexFromItem(items[0]);

  typedef unsigned short MaskPixelType;
  const unsigned int Dimensions = 3;

  typedef itk::Image<MaskPixelType, Dimensions> ImageTypeMask;
  ImageTypeMask::Pointer imageToWrite = convertVtkToItk<MaskPixelType, Dimensions>(mSlicerManagers[index]->mMask);

  typedef itk::MinimumMaximumImageCalculator< ImageTypeMask > CalculatorType;
  CalculatorType::Pointer calculator = CalculatorType::New();
  calculator->SetImage(imageToWrite);
  calculator->Compute();

  if (calculator->GetMaximum() == 0) // this means that the mask image contains no values at all
  {
    ShowErrorMessage("There should be at least one region (near or far) for saving.");
    return;
  }

  auto imageToWrite_wrap = imageToWrite;
  imageToWrite->DisconnectPipeline();
  if (mSlicerManagers[index]->mImageSubType != CAPTK::ImageModalityType::IMAGE_TYPE_PERFUSION)
  {
    ImageTypeMask::DirectionType originalDirection;
    originalDirection[0][0] = mSlicerManagers[index]->mDirection(0, 0);
    originalDirection[0][1] = mSlicerManagers[index]->mDirection(0, 1);
    originalDirection[0][2] = mSlicerManagers[index]->mDirection(0, 2);
    originalDirection[1][0] = mSlicerManagers[index]->mDirection(1, 0);
    originalDirection[1][1] = mSlicerManagers[index]->mDirection(1, 1);
    originalDirection[1][2] = mSlicerManagers[index]->mDirection(1, 2);
    originalDirection[2][0] = mSlicerManagers[index]->mDirection(2, 0);
    originalDirection[2][1] = mSlicerManagers[index]->mDirection(2, 1);
    originalDirection[2][2] = mSlicerManagers[index]->mDirection(2, 2);

    ImageTypeMask::PointType originalOrigin;
    originalOrigin = mSlicerManagers[index]->mOrigin;

    auto infoChanger = itk::ChangeInformationImageFilter< ImageTypeMask >::New();
    infoChanger->SetInput(imageToWrite);
    infoChanger->ChangeDirectionOn();
    infoChanger->ChangeOriginOn();
    infoChanger->SetOutputDirection(originalDirection);
    infoChanger->SetOutputOrigin(originalOrigin);
    infoChanger->Update();
    imageToWrite_wrap = infoChanger->GetOutput();
  }

  QString saveFileName = getSaveFile(this, mInputPathName, mInputPathName + "mask.nii.gz");
  if (!saveFileName.isEmpty())
  {
    std::string filename = saveFileName.toStdString();

    cbica::WriteImage< ImageTypeMask >(imageToWrite_wrap, filename);

    if (cbica::isFile(filename))
    {
      updateProgress(0, "ROI saved!(" + filename + ")");
    }
    else
    {
      ShowErrorMessage("Couldn't write to file: " + filename);
    }
  }
}

void fMainWindow::SaveSeedDrawing()
{
  auto items = m_imagesTable->selectedItems();
  if (items.empty())
    return;

  QString saveFileName = getSaveFile(this, mInputPathName, mInputPathName + "drawing_seed.nii.gz");
  if (!saveFileName.isEmpty())
  {
    std::string filename = saveFileName.toStdString();

    ImageTypeShort3D::Pointer img = convertVtkToItk<ImageTypeShort3D::PixelType, ImageTypeShort3D::ImageDimension>(mSlicerManagers[0]->mMask);
    std::vector<ImageTypeShort3D::IndexType> seedIndices;
    itk::ImageRegionIteratorWithIndex <ImageTypeShort3D> maskIt(img, img->GetLargestPossibleRegion());
    maskIt.GoToBegin();
    while (!maskIt.IsAtEnd())
    {
      if (maskIt.Get() == 3/*seeds are always defined as label '3'*/)
        seedIndices.push_back(maskIt.GetIndex());
      ++maskIt;
    }
    std::string errormsg = "";
    if (seedIndices.size() == 0)
      errormsg = "Draw seed points (label 3) before saving.";
    QMessageBox box(this);
    box.setIcon(QMessageBox::Information);
    box.addButton(QMessageBox::Ok);

    if (errormsg.length() > 0)
    {
      box.setText(QString::fromStdString(errormsg));
      box.setWindowTitle(tr("Error message"));
      box.exec();
      return;
    }

    //save actual near and far region according to the current image type
    std::string InputPixelType = mSlicerManagers[0]->mImage->GetScalarTypeAsString();
    std::string subjectname = m_imagesTable->selectedItems()[0]->data(Qt::DisplayRole).toString().toStdString();
    subjectname = subjectname.substr(0, subjectname.find("_"));
    if (InputPixelType == "short")
    {
      using MaskPixelType = short;
      auto img = convertVtkToItk<MaskPixelType, 3>(mSlicerManagers[0]->mImage);
      mOutputManager.WriteSeedMasks< itk::Image<MaskPixelType, 3> >(img, filename, seedIndices);
    }
    else if (InputPixelType == "unsigned short")
    {
      using MaskPixelType = unsigned short;
      auto img = convertVtkToItk<MaskPixelType, 3>(mSlicerManagers[0]->mImage);
      mOutputManager.WriteSeedMasks<itk::Image<MaskPixelType, 3>>(img, filename, seedIndices);
    }
    else if (InputPixelType == "char")
    {
      using MaskPixelType = char;
      auto img = convertVtkToItk<MaskPixelType, 3>(mSlicerManagers[0]->mImage);
      mOutputManager.WriteSeedMasks<itk::Image<MaskPixelType, 3>>(img, filename, seedIndices);
    }
    else if (InputPixelType == "unsigned char")
    {
      using MaskPixelType = unsigned char;
      auto img = convertVtkToItk<MaskPixelType, 3>(mSlicerManagers[0]->mImage);
      mOutputManager.WriteSeedMasks<itk::Image<MaskPixelType, 3>>(img, filename, seedIndices);
    }
    else if (InputPixelType == "int")
    {
      using MaskPixelType = int;
      auto img = convertVtkToItk<MaskPixelType, 3>(mSlicerManagers[0]->mImage);
      mOutputManager.WriteSeedMasks<itk::Image<MaskPixelType, 3>>(img, filename, seedIndices);
    }
    else if (InputPixelType == "unsigned int")
    {
      using MaskPixelType = unsigned int;
      auto img = convertVtkToItk<MaskPixelType, 3>(mSlicerManagers[0]->mImage);
      mOutputManager.WriteSeedMasks<itk::Image<MaskPixelType, 3>>(img, filename, seedIndices);
    }
    else if (InputPixelType == "double")
    {
      using MaskPixelType = double;
      auto img = convertVtkToItk<MaskPixelType, 3>(mSlicerManagers[0]->mImage);
      mOutputManager.WriteSeedMasks<itk::Image<MaskPixelType, 3>>(img, subjectname, seedIndices);
    }
    else if (InputPixelType == "float")
    {
      using MaskPixelType = float;
      auto img = convertVtkToItk<MaskPixelType, 3>(mSlicerManagers[0]->mImage);
      mOutputManager.WriteSeedMasks<itk::Image<MaskPixelType, 3>>(img, filename, seedIndices);
    }
    else
    {
      cbica::Logging(loggerFile, "Error, input pixel type, '" + InputPixelType + "' is unknown!");
    }

    std::string msg = "Mask saved: " + filename;
    updateProgress(0, msg);

  }
}


void fMainWindow::makeStroke(std::vector<itk::Image<short, 3>::IndexType>& indices, const int value)
{
  std::vector<PointVal> strokePoints;
  for (unsigned int i = 0; i < indices.size(); i++)
  {
    float* pData = (float*)this->mSlicerManagers[0]->GetSlicer(0)->mMask->GetScalarPointer((int)indices[i][0], (int)indices[i][1], (int)indices[i][2]);

    PointVal pt;
    pt.x = indices[i][0];
    pt.y = indices[i][1];
    pt.z = indices[i][1];
    pt.value = (int)*pData;
    strokePoints.push_back(pt);
    *pData = value;
  }
  UpdateAction(strokePoints);
  return;
}
void fMainWindow::clearMask(int label)
{
  if ((mSlicerManagers.size() <= 0))
    return;

  auto img = convertVtkToItk<ImageTypeShort3D::PixelType, 3>(mSlicerManagers[0]->mMask);
  std::vector<ImageTypeShort3D::IndexType> indecesToErase;
  itk::ImageRegionIteratorWithIndex <ImageTypeShort3D> maskIt(img, img->GetLargestPossibleRegion());
  maskIt.GoToBegin();
  while (!maskIt.IsAtEnd())
  {
    if (label < 0)//Clear all
    {
      if (maskIt.Get() > 0)
      {
        indecesToErase.push_back(maskIt.GetIndex());
      }

    }
    else
    {
      if (maskIt.Get() == label)
      {
        indecesToErase.push_back(maskIt.GetIndex());
      }
    }
    ++maskIt;
  }
  makeStroke(indecesToErase, 0);
  this->mSlicerManagers[0]->GetSlicer(0)->mMask->Modified();
  this->mSlicerManagers[0]->Render();
  UpdateNumberOfPointsInTable();
}


//void fMainWindow::StartEGFREstimate()
//{
//}

ImageTypeFloat3D::Pointer fMainWindow::getMaskImage()
{
  ImageTypeFloat3D::Pointer img = NULL;
  if (mSlicerManagers[0]->mMask != NULL)
  {
    img = convertVtkToItk<float, 3>(mSlicerManagers[0]->mMask);
  }

  return img;
}

void fMainWindow::readMaskFile(const std::string &maskFileName)
{
  auto maskFileName_toRead = maskFileName;
  bool imageSanityCheckDone = false;
  if (!mSlicerManagers.empty())
  {
    if (cbica::IsDicom(maskFileName_toRead))
    {
      auto path = cbica::getFilenamePath(maskFileName_toRead);
      auto filesInDir = cbica::filesInDirectory(path, false);

      if (filesInDir.size() == 1) // single DICOM slice
      {
        dicomfilename = maskFileName_toRead;
        maskFileName_toRead = ConversionFrom2Dto3D(maskFileName_toRead);
      }
      else
      {
        auto temp_prev = cbica::normPath(m_tempFolderLocation + "/convertedMask.nii.gz");
        auto maskFromDicom = cbica::ReadImage< ImageTypeFloat3D >(maskFileName_toRead);
        cbica::WriteImage< ImageTypeFloat3D >(maskFromDicom, temp_prev);
        maskFileName_toRead = temp_prev;
      }
    }
    else
    {
      auto maskInfo = cbica::ImageInfo(maskFileName_toRead);
      auto imageSize = mSlicerManagers[0]->mITKImage->GetLargestPossibleRegion().GetSize();
      auto maskSize = maskInfo.GetImageSize();
      if ((imageSize[2] == 1) || (maskSize[2] == 1) || (maskInfo.GetImageDimensions() == 2))
      {
        // this is actually a 2D image which has been loaded in CaPTk as a pseudo-3D image
        auto origin_image = mSlicerManagers[0]->mOrigin;
        auto spacings_image = mSlicerManagers[0]->mITKImage->GetSpacing();
        auto size_image = imageSize;

        auto origin_mask = maskInfo.GetImageOrigins();
        auto spacings_mask = maskInfo.GetImageSpacings();
        auto size_mask = maskInfo.GetImageSize();
        //ImageType::DirectionType directions_image;
        //directions_image[0][0] = mSlicerManagers[0]->mDirection(0, 0);
        //directions_image[0][1] = mSlicerManagers[0]->mDirection(0, 1);
        //directions_image[0][2] = mSlicerManagers[0]->mDirection(0, 2);
        //directions_image[1][0] = mSlicerManagers[0]->mDirection(1, 0);
        //directions_image[1][1] = mSlicerManagers[0]->mDirection(1, 1);
        //directions_image[1][2] = mSlicerManagers[0]->mDirection(1, 2);
        //directions_image[2][0] = mSlicerManagers[0]->mDirection(2, 0);
        //directions_image[2][1] = mSlicerManagers[0]->mDirection(2, 1);
        //directions_image[2][2] = mSlicerManagers[0]->mDirection(2, 2);
        for (size_t i = 0; i < 2; i++)
        {
          if (origin_image[i] != origin_mask[i])
          {
            ShowErrorMessage("The origins of the previously loaded image and mask are inconsistent; cannot load");
            return;
          }
          if (spacings_image[i] != spacings_mask[i])
          {
            auto percentageDifference = std::abs(spacings_image[i] - spacings_mask[i]) * 100;
            percentageDifference /= spacings_image[i];
            if (percentageDifference > 0.0001)
            {
              ShowErrorMessage("The spacings of the previously loaded image and mask are inconsistent; cannot load");
              return;
            }
          }
          if (size_image[i] != size_mask[i])
          {
            ShowErrorMessage("The sizes of the previously loaded image and mask are inconsistent; cannot load");
            return;
          }
        }
        maskFileName_toRead = ConversionFrom2Dto3D(maskFileName_toRead); // all sanity checks passed; load the mask 
        imageSanityCheckDone = true;
      }
    }
    //auto temp_prev = cbica::normPath(m_tempFolderLocation + "/temp_prev.nii.gz");
    auto mask_temp = cbica::ReadImageWithOrientFix< ImageTypeFloat3D >(maskFileName_toRead);
    //SaveImage_withFile(0, temp_prev.c_str());
    if (!imageSanityCheckDone)
    {
      if (!cbica::ImageSanityCheck< ImageTypeFloat3D >(mSlicerManagers[0]->mITKImage, mask_temp))
      {
        ShowErrorMessage("The physical dimensions of the previously loaded image and mask are inconsistent; cannot load");
        return;
      }
      imageSanityCheckDone = true;
    }
    using ImageType = itk::Image<unsigned int, 3>;
    auto inputImage = cbica::ReadImageWithOrientFix< ImageType >(maskFileName_toRead);
    inputImage = ChangeImageDirectionToIdentity< ImageType >(inputImage);

    auto minMaxCalc = itk::MinimumMaximumImageCalculator< ImageType >::New();
    minMaxCalc->SetImage(inputImage);
    minMaxCalc->Compute();
    auto maxVal = minMaxCalc->GetMaximum();

    if (maxVal > 0)
    {
      itk::ImageRegionIteratorWithIndex <ImageType> maskIt(inputImage, inputImage->GetLargestPossibleRegion());
      maskIt.GoToBegin();
      while (!maskIt.IsAtEnd())
      {
        /*
        change to this & also in manual:
        1 for necrosis
        2 for edema
        3 for non-enhancing tumor
        4 for enhancing tumor
        */
        ImageType::IndexType currentIndex = maskIt.GetIndex();
        float* pData = (float*)this->mSlicerManagers[0]->GetSlicer(0)->mMask->GetScalarPointer((int)currentIndex[0], (int)currentIndex[1], (int)currentIndex[2]);
        *pData = 0; // this is done in order to ensure that previously loaded mask is removed
        // this is done to take into account all possible label drawings
        switch (maskIt.Get())
        { // multiLabel: uncomment everything inside this loop and remove references to "near" and "far"
        case DRAW_MODE_LABEL_1:
          *pData = DRAW_MODE_LABEL_1;
          break;
        case 10: // GLISTR defines this as CSF
          *pData = DRAW_MODE_LABEL_7;
          break;
        case DRAW_MODE_LABEL_2:
          *pData = DRAW_MODE_LABEL_2;
          break;
        case 150: // GLISTR defines this is as GM
          *pData = DRAW_MODE_LABEL_5;
          break;
        case DRAW_MODE_LABEL_3:
          *pData = DRAW_MODE_LABEL_3;
          break;
        case 250: // GLISTR defines this is as WM
          *pData = DRAW_MODE_LABEL_3;
          break;
        case DRAW_MODE_LABEL_4:
          *pData = DRAW_MODE_LABEL_4;
          break;
        case 25: // GLISTR defines this is as VS
          *pData = DRAW_MODE_LABEL_4;
          break;
        case DRAW_MODE_LABEL_5: // this is an ambiguous index since GLISTR also uses this for CB
        {
          if (maxVal > DRAW_MODE_LABEL_9) // this means we are reading in GLISTR output
          {
            *pData = DRAW_MODE_LABEL_9;
          }
          else
          {
            *pData = DRAW_MODE_LABEL_5;
          }
          break;
        }
        case 100: // GLISTR defines this is as ED
          *pData = DRAW_MODE_LABEL_2;
          break;
        case DRAW_MODE_LABEL_6:
          *pData = DRAW_MODE_LABEL_6;
          break;
        case 175: // GLISTR defines this is as NCR
          *pData = DRAW_MODE_LABEL_1;
          break;
        case DRAW_MODE_LABEL_7:
          *pData = DRAW_MODE_LABEL_7;
          break;
        case 200: // GLISTR defines this is as TU
          *pData = DRAW_MODE_LABEL_4;
          break;
        case DRAW_MODE_LABEL_8:
          *pData = DRAW_MODE_LABEL_8;
          break;
        case 185: // GLISTR defines this is as NE
          *pData = DRAW_MODE_LABEL_1;
          break;
        case DRAW_MODE_LABEL_9:
          *pData = DRAW_MODE_LABEL_9;
          break;
        case 255: // contingency case in case a map is defined as 255 and 0
          *pData = DRAW_MODE_LABEL_1;
          break;
        default:
          // nothing defined for other cases
          break;
        }
        ++maskIt;
      }
    }
    else
    {
      ShowErrorMessage("Mask file has no pixels greater than '0'");
    }

    UpdateRenderWindows();
    updateProgress(0, "Mask loaded");
  }
  else
  {
    ShowErrorMessage("Please load an image before trying to load an ROI");
    return;
  }

  // Force a render of the mask since updating the render windows doesn't cut it
  this->mSlicerManagers[0]->GetSlicer(0)->mMask->Modified();
  this->mSlicerManagers[0]->Render();
}

std::vector<ImageTypeFloat3D::Pointer> fMainWindow::getLodedImages(std::vector<std::string> &fileNames, std::vector<std::string> &modality, bool onlySelected)
{

  std::vector < ImageTypeFloat3D::Pointer> images;
  if (onlySelected)
  {
    auto items = m_imagesTable->selectedItems();
    if (!items.empty())
    {
      int index = GetSlicerIndexFromItem(items[0]);
      images.push_back(mSlicerManagers[index]->mITKImage);
      fileNames.push_back(mSlicerManagers[index]->mFileName);
      std::string pp = CAPTK::ImageModalityString[mSlicerManagers[index]->mImageSubType];
      modality.push_back(CAPTK::ImageModalityString[mSlicerManagers[index]->mImageSubType]);
    }
  }
  else
  {
    for (unsigned int index = 0; index < mSlicerManagers.size(); index++)
    {
      images.push_back(mSlicerManagers[index]->mITKImage);
      fileNames.push_back(mSlicerManagers[index]->mFileName);
      modality.push_back(CAPTK::ImageModalityString[mSlicerManagers[index]->mImageSubType]);
    }
  }
  return images;
}

void fMainWindow::dragEnterEvent(QDragEnterEvent *event)
{
  if (event->mimeData()->hasFormat("text/uri-list")) {
    event->acceptProposedAction();
  }
}
void fMainWindow::dropEvent(QDropEvent *event)
{
  QList<QUrl> urls = event->mimeData()->urls();
  QStringList vectorOfFiles;
  for (int i = 0; i < (int)urls.size(); i++)
  {
    vectorOfFiles.push_back(urls[i].toLocalFile());
  }
  // if more than 1 files are dropped, assume they are images
  if ((vectorOfFiles.size() > 1) || mSlicerManagers.empty())
  {
    openImages(vectorOfFiles);
  }
  else
  {
    // ask if it is an image or roi
    QMessageBox *box = new QMessageBox(QMessageBox::Question, 
      "Image Type", 
      "Please select the type of image being loaded", QMessageBox::Ok | QMessageBox::Cancel);
    box->button(QMessageBox::Ok)->setText("Image");
    box->button(QMessageBox::Cancel)->setText("ROI");
    box->setAttribute(Qt::WA_DeleteOnClose); //makes sure the msgbox is deleted automatically when closed
    box->setWindowModality(Qt::NonModal);
    QCoreApplication::processEvents();
    if (box->exec() == QMessageBox::Ok)
    {
      openImages(vectorOfFiles);
    }
    else
    {
      readMaskFile(vectorOfFiles[0].toStdString());
    }
  }
}


void fMainWindow::CloseNonViewingDTIImage(QTableWidgetItem* item)
{
  int itemIndexToBeDeleted = 0;
  m_nonVisImagesTable->removeRow(item->row());

  for (unsigned int index = 0; index < mNonViewingImageManager.size(); index++)
  {
    if (mNonViewingImageManager[index]->mImageType == CAPTK::ImageModalityType::IMAGE_TYPE_DTI)
    {
      itemIndexToBeDeleted = index;
      delete mNonViewingImageManager[index];
      break;
    }
  }
  std::vector<SimpleImageManager*>::iterator simpleImageIterator = mNonViewingImageManager.begin();
  for (int i = 0; i < itemIndexToBeDeleted; i++)
    simpleImageIterator++;
  mNonViewingImageManager.erase(simpleImageIterator);
}

void fMainWindow::UpdateNumberOfPointsInTable()
{
  if (mSlicerManagers.size() <= 0)
    return;

  ImageTypeShort3D::Pointer img = convertVtkToItk<short, 3>(mSlicerManagers[0]->mMask);
  int nearCounter = 0;
  int  farCounter = 0;
  int  initCounter = 0;
  typedef itk::ImageRegionIteratorWithIndex <ImageTypeShort3D> IteratorType;
  IteratorType maskIt(img, img->GetLargestPossibleRegion());

  maskIt.GoToBegin();
  while (!maskIt.IsAtEnd())
  {
    if (maskIt.Get() == DRAW_MODE_LABEL_1)
      nearCounter++;
    else if (maskIt.Get() == DRAW_MODE_LABEL_2)
      farCounter++;
    else if (maskIt.Get() == DRAW_MODE_LABEL_3)
      initCounter++;
    ++maskIt;
  }
  mCurrentNearPoints = nearCounter;
  mCurrentFarPoints = farCounter;
  mCurrentInitPoints = initCounter;

}

void fMainWindow::SetPresetComboBox()
{
  presetComboBox->addItem("Auto Scale");
  presetComboBox->addItem("User Scale");
  presetComboBox->addItem("Label Map");
  presetComboBox->addItem("Label Map 2");
  presetComboBox->addItem("Threshold");
  presetComboBox->addItem("Probability");
  presetComboBox->addItem("Geodesic");
}
void fMainWindow::ResetNumberOfPoints()
{
  UpdateNumberOfPointsInTable();
}

ImageTypeFloat3D::Pointer fMainWindow::RescaleImageIntensity(ImageTypeFloat3D::Pointer image)
{
  typedef itk::RescaleIntensityImageFilter< ImageTypeFloat3D, ImageTypeFloat3D > RescaleFilterType;
  RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
  rescaleFilter->SetInput(image);
  rescaleFilter->SetOutputMinimum(0);
  rescaleFilter->SetOutputMaximum(255);
  rescaleFilter->Update();

  ImageTypeFloat3D::Pointer outputimage = rescaleFilter->GetOutput();

  return outputimage;

}
#ifdef BUILD_RECURENCE
void fMainWindow::TrainNewModelOnGivenData(const std::string &inputdirectory, const std::string &outputdirectory, bool useConvData, bool useDTIData, bool usePerfData, bool useDistData)
{
  std::string errorMsg;
  if (inputdirectory.empty())
  {
    ShowErrorMessage("Please provide input directory.");
    help_contextual("Glioblastoma_Recurrence.html");
    return;
  }
  if (outputdirectory.empty())
  {
    ShowErrorMessage("Please provide output directory.");
    help_contextual("Glioblastoma_Recurrence.html");
    return;
  }
  if (!cbica::isDir(outputdirectory))
  {
    if (!cbica::createDir(outputdirectory))
    {
      ShowErrorMessage("Unable to create the output directory");
      help_contextual("Glioblastoma_Recurrence.html");
      return;
    }
  }

  std::vector<double> finalresult;
  std::vector<std::map<CAPTK::ImageModalityType, std::string>> QualifiedSubjects = LoadQualifiedSubjectsFromGivenDirectoryForRecurrence(CAPTK::MachineLearningApplicationSubtype::TRAINING, inputdirectory, useConvData, useDTIData, usePerfData, useDistData);

  if (QualifiedSubjects.size() == 0)
  {
    ShowErrorMessage("The specified directory does not have any subject with all the required imaging sequences.");
    help_contextual("Glioblastoma_Recurrence.html");
    return;
  }

  if (mRecurrenceEstimator.TrainNewModelOnGivenData(QualifiedSubjects, outputdirectory, useConvData, useDTIData, usePerfData, useDistData))
    ShowMessage("Trained infiltration model has been saved at the specified location.", this);
  else
    ShowErrorMessage("Recurrence Estimator wasn't able to save the training files as expected. See log file for details: " + loggerFile);
}
#endif
#ifdef BUILD_PSEUDOPROGRESSION
void fMainWindow::TrainNewPseudoprogressionModelOnGivenData(const std::string &inputdirectory, const std::string &outputdirectory, bool useConvData, bool useDTIData, bool usePerfData, bool useDistData)
{
  std::string errorMsg;
  if (inputdirectory.empty())
  {
    ShowErrorMessage("Please provide input directory.", this);
    help_contextual("Glioblastoma_Pseudoprogression.html");
    return;
  }
  if (outputdirectory.empty())
  {
    ShowErrorMessage("Please provide output directory.", this);
    help_contextual("Glioblastoma_Pseudoprogression.html");
    return;
  }
  if (!cbica::isDir(outputdirectory))
  {
    if (!cbica::createDir(outputdirectory))
    {
      ShowErrorMessage("Unable to create the output directory", this);
      help_contextual("Glioblastoma_Pseudoprogression.html");
      return;
    }
  }

  std::vector<double> finalresult;
  std::vector<std::map<CAPTK::ImageModalityType, std::string>> QualifiedSubjects = LoadQualifiedSubjectsFromGivenDirectoryForPseudoProgression(CAPTK::MachineLearningApplicationSubtype::TRAINING, inputdirectory, useConvData, useDTIData, usePerfData, useDistData);

  if (QualifiedSubjects.size() == 0)
  {
    ShowErrorMessage("The specified directory does not have any subject with all the required imaging sequences.", this);
    help_contextual("Glioblastoma_Pseudoprogression.html");
    return;
  }
  if (QualifiedSubjects.size() > 0 && QualifiedSubjects.size() <= 20)
  {
    ShowErrorMessage("There should be atleast 20 patients to build reliable pseudo-progression model.");
    return;
  }
  if (mPseudoEstimator.TrainNewModelOnGivenData(QualifiedSubjects, outputdirectory, useConvData, useDTIData, usePerfData, useDistData))
    ShowMessage("Trained pseudoprogression model has been saved at the specified location.", this);
  else
    ShowErrorMessage("Pseudoprogression Estimator wasn't able to save the training files as expected. See log file for details: " + loggerFile, this);
}
#endif
#ifdef BUILD_PCA
void fMainWindow::TrainNewPCAModelOnGivenData(const std::string &inputdirectory, const std::string &outputdirectory)
{
  std::string errorMsg;
  if (inputdirectory.empty())
  {
    ShowErrorMessage("Please provide input directory.", this);
    //help_contextual("Glioblastoma_Pseudoprogression.html");
    return;
  }
  if (outputdirectory.empty())
  {
    ShowErrorMessage("Please provide output directory.", this);
    //help_contextual("Glioblastoma_Pseudoprogression.html");
    return;
  }
  if (!cbica::isDir(outputdirectory))
  {
    if (!cbica::createDir(outputdirectory))
    {
      ShowErrorMessage("Unable to create the output directory", this);
      //help_contextual("Glioblastoma_Pseudoprogression.html");
      return;
    }
  }

  std::vector<double> finalresult;
  std::vector<std::map<CAPTK::ImageModalityType, std::string>> QualifiedSubjects = LoadQualifiedSubjectsFromGivenDirectoryForPCA(inputdirectory);

  if (QualifiedSubjects.size() == 0)
  {
    ShowErrorMessage("The specified directory does not have any subject with all the required imaging sequences.", this);
    //help_contextual("Glioblastoma_Pseudoprogression.html");
    return;
  }
  if (QualifiedSubjects.size() > 0 && QualifiedSubjects.size() <= 20)
  {
    ShowErrorMessage("There should be atleast 20 patients to build reliable pseudo-progression model.");
    return;
  }
  PerfusionPCA mPCAEstimator;
  if (mPCAEstimator.PrepareNewPCAModel(10,inputdirectory,outputdirectory,QualifiedSubjects))
    ShowMessage("Trained pseudoprogression model has been saved at the specified location.", this);
  else
    ShowErrorMessage("Pseudoprogression Estimator wasn't able to save the training files as expected. See log file for details: " + loggerFile, this);
}
#endif

//
//void fMainWindow::LoadDicomDrawing()
//{
//  std::string root_directory;
//  QString directory = getExistingDirectory(this, mInputPathName);
//  if (directory.isNull())
//    return;
//
//  typedef itk::Image<unsigned short, 3> InputImageType;
//  typedef itk::ImageSeriesReader< InputImageType > ReaderType;
//  ReaderType::Pointer seriesreader = ReaderType::New();
//
//  typedef itk::GDCMImageIO ImageIOType;
//  ImageIOType::Pointer dicomIO = ImageIOType::New();
//  seriesreader->SetImageIO(dicomIO);
//  typedef itk::GDCMSeriesFileNames NamesGeneratorType;
//  NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();
//  nameGenerator->SetUseSeriesDetails(true);
//  nameGenerator->AddSeriesRestriction("0008|0021");
//
//  nameGenerator->SetInputDirectory(directory.toStdString());
//  try
//  {
//    typedef std::vector< std::string > SeriesIdContainer;
//    const SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();
//
//    SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();
//    SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();
//    while (seriesItr != seriesEnd)
//    {
//      typedef std::vector< std::string > FileNamesContainer;
//      FileNamesContainer fileNames;
//      fileNames = nameGenerator->GetFileNames(seriesItr->c_str());
//      seriesreader->SetFileNames(fileNames);
//      try
//      {
//        seriesreader->Update();
//        typedef unsigned short ROIPixelType;
//        typedef::itk::Image<ROIPixelType, 3> OutputImageType;
//        OutputImageType::Pointer outputImage = seriesreader->GetOutput();
//
//        typedef itk::ImageRegionIteratorWithIndex <OutputImageType> IteratorType;
//        IteratorType maskIt(outputImage, outputImage->GetLargestPossibleRegion());
//        maskIt.GoToBegin();
//
//        auto minMaxCalc = itk::MinimumMaximumImageCalculator< OutputImageType >::New();
//        minMaxCalc->SetImage(outputImage);
//        minMaxCalc->Compute();
//        auto maxVal = minMaxCalc->GetMaximum();
//
//        while (!maskIt.IsAtEnd())
//        {
//          OutputImageType::IndexType currentIndex = maskIt.GetIndex();
//          float* pData = (float*)this->mSlicerManagers[0]->GetSlicer(0)->mMask->GetScalarPointer((int)currentIndex[0], (int)currentIndex[1], (int)currentIndex[2]);
//          // this is done to take into account all possible label drawings
//          switch (maskIt.Get())
//          { // multiLabel: uncomment everything inside this loop and remove references to "near" and "far"
//          case DRAW_MODE_LABEL_1:
//            *pData = DRAW_MODE_LABEL_1;
//            break;
//          case 10: // GLISTR map contingency case
//            *pData = DRAW_MODE_LABEL_1;
//            break;
//          case DRAW_MODE_LABEL_2:
//            *pData = DRAW_MODE_LABEL_2;
//            break;
//          case 150: // GLISTR map contingency case
//            *pData = DRAW_MODE_LABEL_2;
//            break;
//          case DRAW_MODE_LABEL_3:
//            *pData = DRAW_MODE_LABEL_3;
//            break;
//          case 250: // GLISTR map contingency case
//            *pData = DRAW_MODE_LABEL_3;
//            break;
//          case DRAW_MODE_LABEL_4:
//            *pData = DRAW_MODE_LABEL_4;
//            break;
//          case 25: // GLISTR map contingency case
//            *pData = DRAW_MODE_LABEL_4;
//            break;
//          case DRAW_MODE_LABEL_5:
//            if (maxVal > DRAW_MODE_LABEL_9) // if GLISTR map has been defined, this is Cerebellum, i.e., tissue #9
//            {
//              *pData = DRAW_MODE_LABEL_9;
//            }
//            else
//            {
//              *pData = DRAW_MODE_LABEL_5;
//            }
//            break;
//          case 100: // GLISTR map contingency case
//            *pData = DRAW_MODE_LABEL_5;
//            break;
//          case DRAW_MODE_LABEL_6:
//            *pData = DRAW_MODE_LABEL_6;
//            break;
//          case 175: // GLISTR map contingency case
//            *pData = DRAW_MODE_LABEL_6;
//            break;
//          case DRAW_MODE_LABEL_7:
//            *pData = DRAW_MODE_LABEL_7;
//            break;
//          case 200: // GLISTR map contingency case
//            *pData = DRAW_MODE_LABEL_7;
//            break;
//          case DRAW_MODE_LABEL_8:
//            *pData = DRAW_MODE_LABEL_8;
//            break;
//          case 185: // GLISTR map contingency case
//            *pData = DRAW_MODE_LABEL_8;
//            break;
//          case DRAW_MODE_LABEL_9:
//            *pData = DRAW_MODE_LABEL_9;
//            break;
//          default:
//            // nothing defined for other cases
//            break;
//          }
//          ++maskIt;
//        }
//
//        this->mSlicerManagers[0]->GetSlicer(0)->mMask->Modified();
//        this->mSlicerManagers[0]->Render();
//      }
//      catch (itk::ExceptionObject & err)
//      {
//        std::stringstream error;
//        error << err;
//      }
//      ++seriesItr;
//    }
//  }
//  catch (itk::ExceptionObject & err)
//  {
//    std::stringstream error;
//    error << err;
//  }
//}

void fMainWindow::LoadDrawing(const std::string &maskFile)
{
  auto reader = itk::ImageIOFactory::CreateImageIO(maskFile.c_str(), itk::ImageIOFactory::ReadMode);
  if (reader)
  {
    readMaskFile(maskFile);
  }
}

void fMainWindow::LoadDrawing()
{
  if (!mSlicerManagers.empty())
  {
    auto filename = getExistingFile(this, mInputPathName);

    if (filename.isNull() || filename.isEmpty())
    {
      return;
    }
    std::string filename_string = filename.toStdString();
    auto reader = itk::ImageIOFactory::CreateImageIO(filename_string.c_str(), itk::ImageIOFactory::ReadMode);
    if (reader)
    {
      readMaskFile(filename_string);
    }
  }
  else
  {
    ShowErrorMessage("Please load an image before trying to load an ROI", this);
    return;
  }
}

void fMainWindow::UpdateBorderWidget(double startX, double startY, double endX, double endY)
{
  mBorderStartX = std::round(startX);
  mBorderStartY = std::round(startY);
  mBorderEndX = std::round(endX);
  mBorderEndY = std::round(endY);
}
void fMainWindow::UpdateBorderWidget(double startZ, double endZ)
{
  mBorderStartZ = std::round(startZ);
  mBorderStartZ = std::round(startZ);
}

void fMainWindow::overlayUseStateChanged(int state)
{
  if (state == 0)
  {
    for (int i = 0; i < (int)mSlicerManagers.size(); i++)
    {
      for (int j = 0; j < 3; j++) {
        mSlicerManagers[i]->mSlicers[j]->RemoveOverlay();
      }
    }
    UpdateRenderWindows();
  }
}

void fMainWindow::overlaySliderChanged(int value)
{
  auto items = m_imagesTable->selectedItems();
  if (items.empty())
  {
    return;
  }
  int index = GetSlicerIndexFromItem(items[0]);
  if (index >= 0 && index < (int)mSlicerManagers.size())
  {
    for (int i = 0; i < 3; i++)
    {
      mSlicerManagers[index]->GetSlicer(i)->SetOverlayOpacity((double)value / (10 + 1e-6));
    }
  }
  UpdateRenderWindows();
}

void fMainWindow::imageModalityChanged(int value)
{
  for (size_t i = 0; i < mSlicerManagers.size(); i++)
  {
    mSlicerManagers[i]->mImageSubType = imagesPanel->getModality(i);
  }
}

//---------------------------------------------
void fMainWindow::imageSliderChanged()
{
  static int value = -1;
  if (value == image4DSlider->value())
    return;
  else
    value = image4DSlider->value();

  auto items = m_imagesTable->selectedItems();
  if (items.empty())
    return;

  int index = GetSlicerIndexFromItem(items[0]);

  if (mSlicerManagers[index]->mImageSubType == CAPTK::ImageModalityType::IMAGE_TYPE_PERFUSION)
  {
    mSlicerManagers[index]->Get3DImageAtCurrentPerfusionIndex(value);
  }
  AxialViewSliderChanged();
  mSlicerManagers[index]->Picked();
  mSlicerManagers[index]->UpdateViews(0);
  mSlicerManagers[index]->UpdateLinked(0);
  mSlicerManagers[index]->UpdateInfoOnCursorPosition(0);
}
//---------------------------------------------
void fMainWindow::overlayChanged()
{

  overlayChanged(imagesPanel->getSelectedOverlay());
}
void fMainWindow::overlayChanged(QTableWidgetItem *clickedItem)
{
  auto items = m_imagesTable->selectedItems();
  if (items.empty())
  {
    return;
  }
  int slicerManagerIndex = GetSlicerIndexFromItem(items[0]);
  if (slicerManagerIndex < 0 || slicerManagerIndex >= (int)mSlicerManagers.size())
  {
    return;
  }
  //
  int slicerManagerOverlayIndex = GetSlicerIndexFromItem(clickedItem);
  if (slicerManagerOverlayIndex < 0 || slicerManagerOverlayIndex >= (int)mSlicerManagers.size())
  {
    return;
  }
  for (unsigned int i = 0; i < mSlicerManagers.size(); i++)
  {
    if (i != static_cast<unsigned int>(slicerManagerIndex)) {
      for (int j = 0; j < 3; j++) {
        mSlicerManagers[i]->mSlicers[j]->RemoveOverlay();
      }
    }
    else
    {
      for (int j = 0; j < 3; j++)
      {
        mSlicerManagers[slicerManagerIndex]->mSlicers[j]->SetOverlay(mSlicerManagers[slicerManagerOverlayIndex]->mSlicers[j]->mImage);
        //
        double window = mSlicerManagers[slicerManagerOverlayIndex]->mSlicers[j]->GetColorWindow();
        double level = mSlicerManagers[slicerManagerOverlayIndex]->mSlicers[j]->GetColorLevel();
        vtkLookupTable* LUT = static_cast<vtkLookupTable*>(mSlicerManagers[slicerManagerOverlayIndex]->mSlicers[j]->GetWindowLevel()->GetLookupTable());
        if (LUT != NULL)
        {
          mSlicerManagers[slicerManagerIndex]->mSlicers[j]->mOverlayMapper->SetWindow(window);
          mSlicerManagers[slicerManagerIndex]->mSlicers[j]->mOverlayMapper->SetLevel(level);
          mSlicerManagers[slicerManagerIndex]->mSlicers[j]->mOverlayMapper->SetLookupTable(LUT);
        }
        else
        {
          mSlicerManagers[slicerManagerIndex]->mSlicers[j]->mOverlayMapper->SetLookupTable(NULL);
          mSlicerManagers[slicerManagerIndex]->mSlicers[j]->mOverlayMapper->SetWindow(window);
          mSlicerManagers[slicerManagerIndex]->mSlicers[j]->mOverlayMapper->SetLevel(level);
        }
      }
    }
  }
  UpdateRenderWindows();
}

void fMainWindow::openImages(QStringList files, bool callingFromCmd)
{
  if (files.isEmpty())
  {
    if (!callingFromCmd)
    {
      QString extensions = IMAGES_EXTENSIONS;
      extensions += ";;All Files (*)";
      files = QFileDialog::getOpenFileNames(this, tr("Load Images"), mInputPathName, extensions, 0, QFileDialog::DontResolveSymlinks | QFileDialog::DontUseNativeDialog);
      if (files.isEmpty())
        return;
    }
    else
    {
      return;
    }
  }

  int i = 0, fileSizeCheck = files.size() + 1;
  if (mSlicerManagers.empty())
  {
    {
      std::string fileName = files[i].toStdString();
      fileName = cbica::normPath(fileName);
      updateProgress(i + 1, "Opening " + fileName, files.size());
      //auto extension = cbica::getFilenameExtension(fileName);
      //if (!extension.empty())
      //{
      //  std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
      //}
      //if ((extension == ".dcm") || (extension == ".dicom") || (extension == "") ||
      //  (extension == ".ima"))
      if (cbica::IsDicom(fileName))
      {
        QDir d = QFileInfo(fileName.c_str()).absoluteDir();
        QString fname = d.absolutePath();
        dicomfilename = fileName;
        this->openDicomImages(fname);
      }
      else
      {
        LoadSlicerImages(fileName, CAPTK::ImageExtension::NIfTI);
      }
    }
    fileSizeCheck = 1;
  }
  else
  {
    fileSizeCheck = 0;
  }

  // basic sanity check
  if (files.size() > fileSizeCheck)
  {
    std::string erroredFiles, unsupportedExtension;
    std::vector< std::string > basicSanityChecksPassedFiles;
    for (int i = fileSizeCheck; i < files.size(); i++)
    {
      std::string fileName = files[i].toStdString();
      fileName = cbica::normPath(fileName);
      auto extension = cbica::getFilenameExtension(fileName);
      if (!extension.empty())
      {
        std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
      }
      if (!((extension == ".dcm") || (extension == ".dicom") || (extension == "") ||
        (extension == ".ima") || (extension == ".nii") || (extension == ".nii.gz")))
      {
        unsupportedExtension += fileName + "\n";
      }
      else if (!cbica::ImageSanityCheck(files[0].toStdString(), files[i].toStdString()))
      {
        erroredFiles += fileName + "\n";
      }
      else
      {
        basicSanityChecksPassedFiles.push_back(files[i].toStdString());
      }
    }

    if (!unsupportedExtension.empty() && !erroredFiles.empty())
    {
      ShowErrorMessage("Extensions for the following files were not supported, CaPTk will try to load the rest:\n\n" + unsupportedExtension +
        "\n\nAnd the following files are inconsistent with the first loaded image:\n\n" + erroredFiles, this);
      return;
    }
    if (!unsupportedExtension.empty())
    {
      ShowErrorMessage("Extensions for the following files were not supported, CaPTk will try to load the rest:\n\n" + unsupportedExtension, this);
    }
    if (!erroredFiles.empty())
    {
      ShowErrorMessage("Extensions for the following files were not supported, CaPTk will try to load the rest:\n\n" + unsupportedExtension, this);
    }

    for (int i = 0; i < basicSanityChecksPassedFiles.size(); i++)
    {
      std::string fileName = basicSanityChecksPassedFiles[i];
      fileName = cbica::normPath(fileName);
      updateProgress(i + 1, "Opening " + fileName, files.size());
      auto extension = cbica::getFilenameExtension(fileName);
      if (!extension.empty())
      {
        std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
      }
      //if ((extension == ".dcm") || (extension == ".dicom") || (extension == "") ||
      //  (extension == ".ima"))
      if (cbica::IsDicom(fileName))
      {
        QDir d = QFileInfo(fileName.c_str()).absoluteDir();
        QString fname = d.absolutePath();
        dicomfilename = fileName;
        this->openDicomImages(fname);
      }
      else
      {
        LoadSlicerImages(fileName, CAPTK::ImageExtension::NIfTI);
      }
    }
  }

  updateProgress(0, "Loading complete", 100);
}

void fMainWindow::openDicomImages(QString dir)
{
  //QString dir = QFileDialog::getExistingDirectory(this, tr("Open Directory"),
  //  QDir::currentPath(),
  //  QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);

  //if (dir.isNull())
  //{
  //  ShowErrorMessage("Please open a directory containing Dicom images.");
  //  return;
  //}

  //DicomSeriesReader *dicomSeriesReader = new DicomSeriesReader();
  //dicomSeriesReader->SetDirectoryPath(dir.toStdString());
  //bool loadstatus = dicomSeriesReader->LoadDicom();
  //if (!loadstatus)
  //{
  //  QMessageBox::critical(this, "Dicom Loading", "Dicom Load Failed");
  //  return;
  //}

  auto currentImage = cbica::ReadImage<ImageTypeFloat3D>(dir.toStdString());
  if (!currentImage)
  {
    ShowMessage("Dicom Load Failed");
    return;
  }
  SlicerManager* imageManager = new SlicerManager(3, mLandmarks, mSeedPoints, mTissuePoints);
  imageManager->mImageSubType = CAPTK::ImageModalityType::IMAGE_TYPE_UNDEFINED;

  bool bFirstLoad = false;
  if (mSlicerManagers.empty())
  {
    bFirstLoad = true;
  }

  QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));


  imageManager->SetImage(currentImage);
  imageManager->SetOriginalDirection(currentImage->GetDirection());
  imageManager->SetOriginalOrigin(currentImage->GetOrigin());
  //imageManager->SetImage(dicomSeriesReader->GetITKImage());

  //delete dicomSeriesReader; 
  //imageManager->SetFilename(dir.toStdString());
  imageManager->SetMask(mMask);
  imageManager->setTempFolderLocation(m_tempFolderLocation);
  imageManager->mImageSubType = guessImageType(dir.toStdString());
  int rowIndex = (int)mSlicerManagers.size();

  m_imagesTable->setRowCount(rowIndex + 1);
  mSlicerManagers.push_back(imageManager);


  QFileInfo fileinfo(imageManager->GetFileName().c_str());
  std::string seriesDescLabel, seriesDescValue;
  QDir d(dir);
  seriesDescValue = d.dirName().toStdString(); 
  imageManager->SetFilename(seriesDescValue);

  QString id = QString(seriesDescValue.c_str()) + QString::number(mSlicerManagers.size() - 1);
  //
  std::string strImageType = " IMAGE ";

  QTableWidgetItem *item = new QTableWidgetItem(seriesDescValue.c_str());
  item->setData(Qt::UserRole, id.toStdString().c_str());
  item->setFlags(item->flags() & ~Qt::ItemIsEditable);

  QTablePushButton* cButton = new QTablePushButton;
  cButton->setItem(item);
  cButton->setText(QString("X"));
  connect(cButton, SIGNAL(clickedInto(QTableWidgetItem*)), this, SLOT(CloseImage(QTableWidgetItem*)));

  QLabel * label = new QLabel;
  label->setText(QString::fromStdString(strImageType));
  m_imagesTable->setCellWidget(rowIndex, TAB_IMAGES_COLUMN_CLOSE, cButton);
  m_imagesTable->setCellWidget(rowIndex, TAB_IMAGES_COLUMN_TYPE, label);
  m_imagesTable->setItem(rowIndex, TAB_IMAGES_COLUMN_NAME, item);

  imagesPanel->NewImageLoaded(id, imageManager->GetBaseFileName(), rowIndex, strImageType, imageManager->mImageSubType, this);

  mSlicerManagers.back()->SetId(id.toStdString());
  connect(mSlicerManagers.back(), SIGNAL(LeftButtonReleaseSignal(int)), this, SLOT(propogateSlicerPosition(int)));
  connect(mSlicerManagers.back(), SIGNAL(currentImageChanged(std::string &)), this, SLOT(CurrentImageChanged(std::string &)));
  connect(mSlicerManagers.back(), SIGNAL(currentPickedImageChanged(std::string)), this, SLOT(CurrentPickedImageChanged(std::string)));
  connect(mSlicerManagers.back(), SIGNAL(UpdatePosition(int, double, double, double, double, double, double, double)), this, SLOT(MousePositionChanged(int, double, double, double, double, double, double, double)));
  connect(mSlicerManagers.back(), SIGNAL(WindowLevelChanged()), this, SLOT(WindowLevelChanged()));
  connect(mSlicerManagers.back(), SIGNAL(UpdateSlice(int, int)), this, SLOT(UpdateSlice(int, int)));
  connect(mSlicerManagers.back(), SIGNAL(UpdateSliceRange(int, int, int)), this, SLOT(UpdateSliceRange(int, int, int)));
  connect(mSlicerManagers.back(), SIGNAL(UpdateLinkManager(std::string, int, double, double, double)), this, SLOT(UpdateLinkManager(std::string, int, double, double, double)));
  connect(mSlicerManagers.back(), SIGNAL(ChangeImageWithOrder(SlicerManager*, int)), this, SLOT(ChangeImageWithOrder(SlicerManager*, int)));
  connect(mSlicerManagers.back(), SIGNAL(UpdateBorderWidgetInMain(double, double, double, double)), this, SLOT(UpdateBorderWidget(double, double, double, double)));
  connect(mSlicerManagers.back(), SIGNAL(UpdateBorderWidgetInMain(double, double)), this, SLOT(UpdateBorderWidget(double, double)));
  connect(mSlicerManagers.back(), SIGNAL(UpdateActionInMain(const QVariantList&)), this, SLOT(UpdateActionQ(const QVariantList&)));

  //connect(mSlicerManagers.back(), SIGNAL(SeedPointsAdded()), tumorPanel, SLOT(sAddPoint()));
  //connect(mSlicerManagers.back(), SIGNAL(SeedPointsAdded(int, bool)), tumorPanel, SLOT(sAddPoint(int, bool)));
  //connect(mSlicerManagers.back(), SIGNAL(TissuePointsAdded(int)), tumorPanel, SLOT(tAddPoint(int)));
  //connect(m_tabWidget, SIGNAL(currentChanged(int)), tumorPanel, SLOT(tabSelected()));
  InitSlicers();

  if (bFirstLoad)
  {
    InitMask(mSlicerManagers.back()->mImage);
  }
  for (int j = 0; j < (int)mSlicerManagers.back()->mSlicers.size(); j++)
  {
    mSlicerManagers.back()->mSlicers[j]->SetMask(mSlicerManagers.back()->GetMask());
  }

  if (mSlicerManagers.size() > 0)
  {
    if (mSlicerManagers.back()->mMask->GetDimensions()[2] != 1)
    {
      CoronalViewWidget->show();
      SaggitalViewWidget->show();
    }
    AxialViewWidget->show();
    infoPanel->show();

    windowLabel->setEnabled(true);
    windowSpinBox->setEnabled(true);
    levelLabel->setEnabled(true);
    levelSpinBox->setEnabled(true);
    presetLabel->setEnabled(true);
    presetComboBox->setEnabled(true);

    if (bFirstLoad)
    {
      for (int i = 0; i < 3; i++)
      {
        mSlicerManagers.back()->GetSlicer(i)->SetInitPosition();
      }
      QTableWidgetItem* item = NULL;
      item = GetItemFromSlicerManager(mSlicerManagers[0]);
      DisplayChanged(item);
    }
    else
    {
      QTableWidgetItem* item = NULL;
      for (int i = 0; i < (int)mSlicerManagers.size(); i++)
      {
        item = GetItemFromSlicerManager(mSlicerManagers[i]);
        if (!item->isSelected())
        {
          item->setSelected(true);
        }
      }
      DisplayChanged(item);
    }

    if (mSlicerManagers.size() > 1)
    {
      for (int i = 0; i < (int)mSlicerManagers.size(); i++)
      {
        for (int j = i + 1; j < (int)mSlicerManagers.size(); j++)
        {
          AddLink(/*QString::fromStdString*/(mSlicerManagers[i]->GetId().c_str()), /*QString::fromStdString*/(mSlicerManagers[j]->GetId().c_str()));
        }
      }
    }
    QTableWidgetItem* item = GetItemFromSlicerManager(mSlicerManagers.back());
    item->setSelected(true);
    InitDisplay();
  }
  propogateSlicerPosition();
  updateProgress(0);
  QApplication::restoreOverrideCursor();

}

void fMainWindow::ApplicationGeodesicTreshold()
{
  if (m_imgGeodesicOut.IsNull())
  {
    return;
  }
  itk::ImageRegionIterator<ImageTypeShort3D> imageIterator(m_imgGeodesicOut, m_imgGeodesicOut->GetLargestPossibleRegion());
  while (!imageIterator.IsAtEnd())
  {
    auto currentIndex = imageIterator.GetIndex();
    float val = imageIterator.Get();
    float* pData = (float*)this->mSlicerManagers[0]->GetSlicer(0)->mMask->GetScalarPointer((int)currentIndex[0], (int)currentIndex[1], (int)currentIndex[2]);
    if (val < thresholdSpinBox->value())
    {
      *pData = 1.0;
    }
    else
    {
      *pData = 0.0;
    }
    ++imageIterator;
  }
  this->mSlicerManagers[0]->GetSlicer(0)->mMask->Modified();
  this->mSlicerManagers[0]->Render();
}
void fMainWindow::ImageDenoising()
{
  auto items = m_imagesTable->selectedItems();
  if (items.empty())
  {
    ShowErrorMessage("Please load the image you would like to de-noise", this);
    return;
  }

  int index = GetSlicerIndexFromItem(items[0]);
  if (index < 0 || index >= (int)mSlicerManagers.size())
    return;

  QString saveFileName = getSaveFile(this, mInputPathName, mInputPathName + "denoise.nii.gz");
  if (!saveFileName.isEmpty())
  {
    auto saveFileName_string = saveFileName.toStdString();
    //Job TBD replace with app name
    typedef itk::ImageDuplicator<ImageTypeFloat3D> DuplicatorType;
    DuplicatorType::Pointer duplicator = DuplicatorType::New();
    duplicator->SetInputImage(mSlicerManagers[index]->mITKImage);
    duplicator->Update();
    ImageTypeFloat3D::Pointer inputImage = duplicator->GetOutput();
    updateProgress(5, "Susan noise removal in process");
    SusanDenoising denoising /*= SusanDenoising()*/;
    ImageTypeFloat3D::Pointer outputImage = denoising.Run<ImageTypeFloat3D>(inputImage);
    if (outputImage.IsNotNull())
    {
      updateProgress(80, "Saving file");
      cbica::WriteImage< ImageTypeFloat3D >(outputImage, saveFileName_string);
      if (cbica::fileExists(saveFileName_string))
      {
        updateProgress(90, "Displaying output");
        LoadSlicerImages(saveFileName_string, CAPTK::ImageExtension::NIfTI);

      }
      updateProgress(0, "Susan noise removal finished");
    }
    else
    {
      updateProgress(0, "Error in Susan noise removal!!");
    }
  }
}


void fMainWindow::ImageBiasCorrection()
{
  auto items = m_imagesTable->selectedItems();
  if (items.empty())
  {
    ShowErrorMessage("Please load an image to run bias correction on", this);
    return;
  }

  int index = GetSlicerIndexFromItem(items[0]);
  if (index < 0 || index >= (int)mSlicerManagers.size())
    return;

  QString saveFileName = getSaveFile(this, mInputPathName, mInputPathName + "biasCorrect.nii.gz");
  if (!saveFileName.isEmpty())
  {
    auto saveFileName_string = saveFileName.toStdString();
    typedef itk::ImageDuplicator <ImageTypeFloat3D > DuplicatorType;
    DuplicatorType::Pointer duplicator = DuplicatorType::New();
    duplicator->SetInputImage(mSlicerManagers[index]->mITKImage);
    duplicator->Update();
    auto inputImage = duplicator->GetOutput();
    updateProgress(5, "Bias correction in process");
    N3BiasCorrection biasCorrecter /*= N3BiasCorrection()*/;
    auto outputImage = biasCorrecter.Run<ImageTypeFloat3D>(inputImage);
    if (outputImage.IsNotNull())
    {
      updateProgress(80, "Saving file");
      cbica::WriteImage< ImageTypeFloat3D >(outputImage, saveFileName_string);
      if (cbica::fileExists(saveFileName_string))
      {
        updateProgress(90, "Displaying output");
        LoadSlicerImages(saveFileName_string, CAPTK::ImageExtension::NIfTI);
      }
      updateProgress(0, "Bias correction finished");
    }
    else
    {
      updateProgress(0, "Error in Bias correction!");
    }
  }
}
void fMainWindow::ImageRegistration()
{
  registrationPanel.mInputPathName = mInputPathName;
  registrationPanel.exec();
}

void fMainWindow::ImageHistogramMatching()
{
  // open a simple dialog box with reference image, input and output
  histoMatchPanel.SetCurrentImagePath(mInputPathName);
  histoMatchPanel.exec();
}

void fMainWindow::ImageDeepMedicNormalizer()
{
#ifndef WIN32
  ShowErrorMessage("DeepMedic is currently not available for your platform but will be soon.", this);
  return;
#endif
  // open a simple dialog box with reference image, input and output
  //deepMedicNormPanel.exec();
}

void fMainWindow::ImageSkullStripping()
{
  // open a simple dialog box with reference image, input and output
  //skullStrippingPanel.exec();
}

void fMainWindow::ApplicationTheia()
{
  if (!mSlicerManagers.empty())
  {
    if (isMaskDefined())
    {
      std::string maskFile = m_tempFolderLocation + "/theia_mask.nii.gz";
      cbica::WriteImage< ImageTypeFloat3D >(getMaskImage(), maskFile);

      auto items = m_imagesTable->selectedItems();
      auto index = GetSlicerIndexFromItem(items[0]);

      QStringList args;
      args << "-i" << mSlicerManagers[index]->GetFileName().c_str() << "-m" << maskFile.c_str();
      startExternalProcess(getApplicationPath("Theia").c_str(), args);
    }
    else
    {
      ShowErrorMessage("Please initialize a valid mask before trying 3D Visualizer", this);
      return;
    }
  }
  else
  {
    ShowErrorMessage("Please load at least a single image before trying 3D Visualizer", this);
    return;
  }
}

void fMainWindow::EnableComparisonMode(bool enable)
{
  int nLoadedData = mSlicerManagers.size();
  if (nLoadedData < 2 || nLoadedData > 3)
  {
    ShowMessage("Comparison mode only works with 2 or 3 datasets. Please load 2 or 3 datasets to enable comparison mode", this);
    return;
  }
  if ((mSlicerManagers[0]->mITKImage->GetLargestPossibleRegion().GetSize()[2] == 1)) //! e.g. Mammography images
  {
    ShowErrorMessage("2D images are not currently supported in Comparison Mode.");
    return;
  }

  this->SetComparisonMode(enable);

  if (enable) //! enabling comparison
  {
    if (m_ComparisonViewerLeft.GetPointer() == nullptr &&
      m_ComparisonViewerCenter.GetPointer() == nullptr &&
      m_ComparisonViewerRight.GetPointer() == nullptr)
    {
      m_ComparisonViewerLeft = vtkSmartPointer<Slicer>::New();
      m_ComparisonViewerCenter = vtkSmartPointer<Slicer>::New();
      m_ComparisonViewerRight = vtkSmartPointer<Slicer>::New();

	  for (int i = 0; i < this->GetComparisonViewers().size(); i++)
	  {
		  this->GetComparisonViewers()[i]->SetComparisonMode(true);
	  }
    }

	for (int i = 0; i < this->GetComparisonViewers().size(); i++)
	{
		this->GetComparisonViewers()[i]->SetImage(mSlicerManagers[i]->GetSlicer(0)->GetImage(), mSlicerManagers[i]->GetSlicer(0)->GetTransform());
		this->GetComparisonViewers()[i]->SetMask(mSlicerManagers[0]->GetMask());
		this->GetComparisonViewers()[i]->SetRenderWindow(0, nullptr);
		this->GetComparisonViewers()[i]->SetImageSeriesDescription(mSlicerManagers[i]->mBaseFileName);
	}

	if (nLoadedData == 2) //! 2 datasets are loaded
	{
		m_ComparisonViewerLeft->SetRenderWindow(0, AxialViewWidget->GetRenderWindow());
		m_ComparisonViewerCenter->SetRenderWindow(0, CoronalViewWidget->GetRenderWindow());

		SaggitalViewWidget->hide();
		SaggitalViewSlider->hide();
	}
	else if (nLoadedData == 3) //! 3 datasets are loaded
	{
		m_ComparisonViewerLeft->SetRenderWindow(0, AxialViewWidget->GetRenderWindow());
		m_ComparisonViewerCenter->SetRenderWindow(0, CoronalViewWidget->GetRenderWindow());
		m_ComparisonViewerRight->SetRenderWindow(0, SaggitalViewWidget->GetRenderWindow());
	}

      for (int i = 0; i < this->GetComparisonViewers().size(); i++)
      {
        InteractorStyleNavigator* style = InteractorStyleNavigator::New();
        ComparisonViewerCommand *smc = ComparisonViewerCommand::New();
        smc->SetCurrentViewer(this->GetComparisonViewers()[i]);
        smc->SetComparisonViewers(this->GetComparisonViewers());
        smc->SM = mSlicerManagers[0];
        style->AddObserver(vtkCommand::KeyPressEvent, smc);
        style->AddObserver(vtkCommand::WindowLevelEvent, smc);
        style->AddObserver(vtkCommand::EndWindowLevelEvent, smc);
        style->AddObserver(vtkCommand::StartWindowLevelEvent, smc);
        style->AddObserver(vtkCommand::PickEvent, smc);
        style->AddObserver(vtkCommand::StartPickEvent, smc);
        style->AddObserver(vtkCommand::LeaveEvent, smc);
        style->AddObserver(vtkCommand::UserEvent, smc);
        style->AddObserver(vtkCommand::MouseWheelForwardEvent, smc);
        style->AddObserver(vtkCommand::MouseWheelBackwardEvent, smc);
        style->AddObserver(vtkCommand::LeftButtonReleaseEvent, smc);
        style->AddObserver(vtkCommand::EndPickEvent, smc);
        style->AddObserver(vtkCommand::EndInteractionEvent, smc);
        style->SetAutoAdjustCameraClippingRange(1);
        this->GetComparisonViewers()[i]->SetInteractorStyle(style);
        style->Delete();
      }

      //! when we enter comparison mode, the WL should be same as in regular mode
      std::vector<vtkSmartPointer<Slicer>> comparisonViewers = this->GetComparisonViewers();
      for (int i = 0; i < comparisonViewers.size(); i++)
      {
        comparisonViewers[i]->SetColorWindow(windowSpinBox->value());
        comparisonViewers[i]->SetColorLevel(levelSpinBox->value());
      }

	  for (int i = 0; i < comparisonViewers.size(); i++)
	  {
		  comparisonViewers[i]->SetDisplayMode(true);
	  }

      //!comparison mode connections
      disconnect(AxialViewSlider, SIGNAL(valueChanged(int)), this, SLOT(AxialViewSliderChanged()));
      disconnect(CoronalViewSlider, SIGNAL(valueChanged(int)), this, SLOT(CoronalViewSliderChanged()));
      disconnect(SaggitalViewSlider, SIGNAL(valueChanged(int)), this, SLOT(SaggitalViewSliderChanged()));

      connect(AxialViewSlider, SIGNAL(valueChanged(int)), this, SLOT(OnSliderMovedInComparisonMode(int)));
      connect(CoronalViewSlider, SIGNAL(valueChanged(int)), this, SLOT(OnSliderMovedInComparisonMode(int)));
      connect(SaggitalViewSlider, SIGNAL(valueChanged(int)), this, SLOT(OnSliderMovedInComparisonMode(int)));

	  for (int i = 0; i < comparisonViewers.size(); i++)
	  {
		  comparisonViewers[i]->Render();
	  }
  }
  else
  {
	  //! disabling comparison and coming back to regular mode

	  if (nLoadedData == 2) //! 2 datasets loaded
	  {
		  mSlicerManagers[0]->SetImage(mSlicerManagers[0]->GetITKImage());
		  mSlicerManagers[1]->SetImage(mSlicerManagers[1]->GetITKImage());

		  mSlicerManagers[0]->GetSlicer(0)->SetRenderWindow(0, nullptr);
		  mSlicerManagers[1]->GetSlicer(0)->SetRenderWindow(0, nullptr);

		  mSlicerManagers[0]->GetSlicer(0)->SetRenderWindow(0, AxialViewWidget->GetRenderWindow());
		  mSlicerManagers[1]->GetSlicer(0)->SetRenderWindow(0, AxialViewWidget->GetRenderWindow());

		  SaggitalViewWidget->show();
		  SaggitalViewSlider->show();
	  }
	  else if (nLoadedData == 3) //! 3 datasets loaded
	  {
		  mSlicerManagers[0]->SetImage(mSlicerManagers[0]->GetITKImage());
		  mSlicerManagers[1]->SetImage(mSlicerManagers[1]->GetITKImage());
		  mSlicerManagers[2]->SetImage(mSlicerManagers[2]->GetITKImage());

		  mSlicerManagers[0]->GetSlicer(0)->SetRenderWindow(0, nullptr);
		  mSlicerManagers[1]->GetSlicer(0)->SetRenderWindow(0, nullptr);
		  mSlicerManagers[2]->GetSlicer(0)->SetRenderWindow(0, nullptr);

		  mSlicerManagers[0]->GetSlicer(0)->SetRenderWindow(0, AxialViewWidget->GetRenderWindow());
		  mSlicerManagers[1]->GetSlicer(0)->SetRenderWindow(0, AxialViewWidget->GetRenderWindow());
		  mSlicerManagers[2]->GetSlicer(0)->SetRenderWindow(0, AxialViewWidget->GetRenderWindow());
	  }

    //!regular mode connections
    connect(AxialViewSlider, SIGNAL(valueChanged(int)), this, SLOT(AxialViewSliderChanged()));
    connect(CoronalViewSlider, SIGNAL(valueChanged(int)), this, SLOT(CoronalViewSliderChanged()));
    connect(SaggitalViewSlider, SIGNAL(valueChanged(int)), this, SLOT(SaggitalViewSliderChanged()));

    disconnect(AxialViewSlider, SIGNAL(valueChanged(int)), this, SLOT(OnSliderMovedInComparisonMode(int)));
    disconnect(CoronalViewSlider, SIGNAL(valueChanged(int)), this, SLOT(OnSliderMovedInComparisonMode(int)));
    disconnect(SaggitalViewSlider, SIGNAL(valueChanged(int)), this, SLOT(OnSliderMovedInComparisonMode(int)));

	for (int i = 0; i < this->GetComparisonViewers().size(); i++)
	{
		this->GetComparisonViewers()[i]->SetDisplayMode(false);
	}

    this->InitDisplay();

    mSlicerManagers[0]->Render();

  }
}

void fMainWindow::ApplicationDeepMedicSegmentation(int type)
{
  if (type <= fDeepMedicDialog::SkullStripping) // different cases for individual models can be put in this way
  {
    if (mSlicerManagers.size() < 4)
    {
      ShowErrorMessage("This model needs the following images to work: T1CE, T1, T2, FLAIR", this);
      return;
    }
  }

  // redundancy check
  if (type >= fDeepMedicDialog::Max)
  {
    ShowErrorMessage("Unsupported model type, please check", this);
    return;
  }

  deepMedicDialog.SetDefaultModel(type);
  deepMedicDialog.SetCurrentImagePath(mInputPathName);
  deepMedicDialog.exec();
}

void fMainWindow::CallDeepMedicSegmentation(const std::string modelDirectory, const std::string outputDirectory)
{
  std::string file_t1ce, file_t1, file_flair, file_t2;

  cbica::createDir(outputDirectory);

  auto file_mask = outputDirectory + "/dm_mask.nii.gz";
  if (!isMaskDefined())
  {
    file_mask = "";
  }
  else
  {
    cbica::WriteImage< TImageType >(getMaskImage(), file_mask);
  }

  if (!cbica::isDir(modelDirectory))
  {
    ShowErrorMessage("Model directory was not found, please try with another");
    return;
  }
  if (!cbica::isFile(modelDirectory + "/modelConfig.txt"))
  {
    ShowErrorMessage("'modelConfig.txt' was not found in the directory, please check");
    return;
  }
  //if (!cbica::isFile(modelDirectory + "/model.ckpt"))
  //{
  //  ShowErrorMessage("'model.ckpt' was not found in the directory, please check");
  //  return;
  //}

  auto modelConfigFile = modelDirectory + "/modelConfig.txt",
    modelCkptFile = modelDirectory + "/model.ckpt";

  std::string files_forCommand;

  int progressBar = 0;
  for (size_t i = 0; i < mSlicerManagers.size(); i++)
  {
    switch (mSlicerManagers[i]->mImageSubType)
    {
    case CAPTK::ImageModalityType::IMAGE_TYPE_T1CE:
    {
      auto temp = cbica::normPath(m_tempFolderLocation + "/t1ce.nii.gz");
      SaveImage_withFile(i, temp.c_str());
      file_t1ce = temp;
      break;
    }
    case CAPTK::ImageModalityType::IMAGE_TYPE_T1:
    {
      auto temp = cbica::normPath(m_tempFolderLocation + "/t1.nii.gz");
      SaveImage_withFile(i, temp.c_str());
      file_t1 = temp;
      break;
    }
    case CAPTK::ImageModalityType::IMAGE_TYPE_T2:
    {
      auto temp = cbica::normPath(m_tempFolderLocation + "/t2.nii.gz");
      SaveImage_withFile(i, temp.c_str());
      file_t2 = temp;
      break;
    }
    case CAPTK::ImageModalityType::IMAGE_TYPE_T2FLAIR:
    {
      auto temp = cbica::normPath(m_tempFolderLocation + "/flair.nii.gz");
      SaveImage_withFile(i, temp.c_str());
      file_flair = temp;
      break;
    }
    default:
      ShowErrorMessage("DeepMedic needs the following images to work: T1-Gd, T1, T2, FLAIR", this);
      break;
    }
  }

  QMessageBox *box = new QMessageBox(QMessageBox::Question, "Long running Application",
    "Deep Learning inference takes 5-30 minutes to run, during which FeTS will not be responsive; press OK to continue...",
    QMessageBox::Ok | QMessageBox::Cancel);
  box->setAttribute(Qt::WA_DeleteOnClose); //makes sure the msgbox is deleted automatically when closed
  box->setWindowModality(Qt::NonModal);
  QCoreApplication::processEvents();
  if (box->exec() == QMessageBox::Ok)
  {
    // TBD: this requires cleanup
    int type;
    if (modelDirectory.find("tumor") != std::string::npos)
    {
      type = 0;
    }
    else if (modelDirectory.find("skull") != std::string::npos)
    {
      type = 1;
    }

    QStringList args;
    args << "-md" << modelDirectory.c_str() << "-o" << outputDirectory.c_str();

    // parsing the modality-agnostic case
    auto modelDir_lower = modelDirectory;
    std::transform(modelDir_lower.begin(), modelDir_lower.end(), modelDir_lower.begin(), ::tolower);
    if (modelDir_lower.find("modalityagnostic") != std::string::npos)
    {
      // we only want to pick up a single modality, in this case, so the first one loaded is picked
      // order of preference is t1, t1ce, t2, fl
      if (!file_t1.empty())
      {
        files_forCommand += file_t1 + ",";
      }
      else if (!file_t1ce.empty())
      {
        files_forCommand += file_t1ce + ",";
      }
      else if (!file_t2.empty())
      {
        files_forCommand += file_t2 + ",";
      }
      else if (!file_flair.empty())
      {
        files_forCommand += file_flair + ",";
      }
    }
    else
    {
      files_forCommand = file_t1 + "," + file_t1ce + "," + file_t2 + "," + file_flair + ",";
    }
    files_forCommand.pop_back(); // last "," removed

    args << "-i" << files_forCommand.c_str() << "-o" << outputDirectory.c_str();

    if (!file_mask.empty())
    {
      args << "-m" << file_mask.c_str();
    }
    updateProgress(5, "Starting DeepMedic Segmentation");

    auto dmExe = getApplicationPath("DeepMedic");
    if (!cbica::exists(dmExe))
    {
      ShowErrorMessage("DeepMedic executable doesn't exist; can't run");
      updateProgress(0, "");
      return;
    }


    if (startExternalProcess(dmExe.c_str(), args) != 0)
    {
      ShowErrorMessage("DeepMedic returned with exit code != 0");
      updateProgress(0, "");
      return;
    }

    auto output = outputDirectory + "/predictions/testApiSession/predictions/Segm.nii.gz";
    if (cbica::exists(output))
    {
      readMaskFile(output);
      updateProgress(100, "Completed.");
    }
    else
    {
      ShowErrorMessage("DeepMedic failed to generate results");
      updateProgress(0, "");
    }
  }

  return;
}

void fMainWindow::CallImageSkullStripping(const std::string referenceAtlas, const std::string referenceMask,
  const std::string inputImageFile, const std::string outputImageFile)
{
  ShowErrorMessage("Skull Stripping takes a long time to run, during which CaPTk will not be responsive.", this, "Long Running Application");
  if (!cbica::isFile(referenceAtlas))
  {
    ShowErrorMessage("Reference Atlas is not a valid file, please re-check", this);
    return;
  }
  if (!cbica::isFile(referenceMask))
  {
    ShowErrorMessage("Reference Mask is not a valid file, please re-check", this);
    return;
  }
  if (!cbica::isFile(inputImageFile))
  {
    ShowErrorMessage("Input Image is not a valid file, please re-check", this);
    return;
  }
  auto referenceAtlasImage = cbica::ReadImage< ImageTypeFloat3D >(referenceAtlas);
  auto referenceAtlasMaskImage = cbica::ReadImage< ImageTypeFloat3D >(referenceMask);
  auto inputImageImage = cbica::ReadImage< ImageTypeFloat3D >(inputImageFile);

  auto outputImage = cbica::GetSkullStrippedImage< ImageTypeFloat3D >(inputImageImage, referenceAtlasImage, referenceAtlasMaskImage);

  if ((cbica::getFilenameExtension(outputImageFile) != ".nii") && (cbica::getFilenameExtension(outputImageFile) != ".nii.gz"))
  {
    std::string path, base, ext;
    cbica::splitFileName(outputImageFile, path, base, ext);
    cbica::WriteImage< ImageTypeFloat3D >(outputImage, path + base + ".nii.gz");
    LoadSlicerImages(path + base + ".nii.gz", CAPTK::ImageExtension::NIfTI);
  }
  else
  {
    cbica::WriteImage< ImageTypeFloat3D >(outputImage, outputImageFile);
    LoadSlicerImages(outputImageFile, CAPTK::ImageExtension::NIfTI);
  }
}

void fMainWindow::CallLabelValuesChange(const std::string oldValues, const std::string newValues)
{
  if (!isMaskDefined())
  {
    ShowErrorMessage("A valid mask needs to be loaded");
    return;
  }
  auto oldValues_string_split = cbica::stringSplit(oldValues, "x");
  auto newValues_string_split = cbica::stringSplit(newValues, "x");

  if (oldValues_string_split.size() != newValues_string_split.size())
  {
    ShowErrorMessage("Old and New values have the same number of inputs", this);
    return;
  }

  auto output = cbica::ChangeImageValues< ImageTypeFloat3D >(getMaskImage(), oldValues, newValues);

  if (output.IsNull())
  {
    ShowErrorMessage("Changing values did not work as expected, please try again with correct syntax");
    return;
  }

  std::string tempFile = m_tempFolderLocation + "/mask_changedValues.nii.gz";
  cbica::WriteImage< ImageTypeFloat3D >(output, tempFile);
  readMaskFile(tempFile);
}

void fMainWindow::CallImageHistogramMatching(const std::string referenceImage, const std::string inputImageFile, const std::string outputImageFile)
{
  if (!referenceImage.empty() && !inputImageFile.empty() && !outputImageFile.empty())
  {
    if (!cbica::isFile(referenceImage))
    {
      ShowErrorMessage("Reference Image is not a valid file, please re-check", this);
      return;
    }
    if (!cbica::isFile(inputImageFile))
    {
      ShowErrorMessage("Input Image is not a valid file, please re-check", this);
      return;
    }
    auto referenceAtlasImage = cbica::ReadImage< ImageTypeFloat3D >(referenceImage);
    auto inputImageImage = cbica::ReadImage< ImageTypeFloat3D >(inputImageFile);

    auto outputImage = cbica::GetHistogramMatchedImage< ImageTypeFloat3D >(inputImageImage, referenceAtlasImage);

    cbica::WriteImage< ImageTypeFloat3D >(outputImage, outputImageFile);

    LoadSlicerImages(outputImageFile, CAPTK::ImageExtension::NIfTI);
  }
  else
  {
    ShowErrorMessage("Please provide all inputs before trying histogram matching", this);
    help_contextual("preprocessing_histoMatch.html");
    return;
  }
}

void fMainWindow::CustomPreprocessing()
{

}
void fMainWindow::ChangeBrushSize(int size)
{
  updateDrawMode();
}

void fMainWindow::ChangeMaskOpacity(int newMaskOpacity) // multiLabel uncomment this function
{
  double tempOpacity = newMaskOpacity * 0.1;
  for (size_t i = 0; i < this->mSlicerManagers.size(); i++)
  {
    for (size_t j = 0; j < 3; j++)
    {
      this->mSlicerManagers[i]->GetSlicer(j)->mMaskOpacity = tempOpacity;
      this->mSlicerManagers[i]->GetSlicer(j)->mMaskActor->SetOpacity(tempOpacity);
      this->mSlicerManagers[i]->GetSlicer(j)->mMask->Modified();
    }
  }
  this->mSlicerManagers[0]->Render();
}

void fMainWindow::ChangeDrawingLabel(int drawingLabel) // multiLabel uncomment this function
{
  updateDrawMode();
}

/**
* Read a transform specification, format file,number
*/
TransformSpec read_transform_spec(std::string &file)
{
  std::string spec = file;
  size_t pos = spec.find_first_of(',');

  TransformSpec ts;
  ts.filename = spec.substr(0, pos);
  ts.exponent = 1.0;

  if (!itksys::SystemTools::FileExists(ts.filename.c_str()))
    throw GreedyException("File '%s' does not exist", ts.filename.c_str());

  if (pos != std::string::npos)
  {
    errno = 0; char *pend;
    std::string expstr = spec.substr(pos + 1);
    ts.exponent = std::strtod(expstr.c_str(), &pend);

    if (errno || *pend)
      throw GreedyException("Expected a floating point number after comma in transform specification, instead got '%s'",
        spec.substr(pos).c_str());

  }
  return ts;
}

//Reads radius for registration
std::vector<int> read_int_vector(std::string &nccRadii)
{
  std::string arg = nccRadii;
  std::istringstream f(arg);
  std::string s;
  std::vector<int> vector;
  while (getline(f, s, 'x'))
  {
    errno = 0; char *pend;
    long val = std::strtol(s.c_str(), &pend, 10);
    //std::cout << "Radii: " << val <<std::endl;
    if (errno || *pend)
      throw GreedyException("Expected an integer vector as parameter, instead got '%s'",
        arg.c_str());
    vector.push_back((int)val);
  }

  if (!vector.size())
    throw GreedyException("Expected an integer vector as parameter, instead got '%s'",
      arg.c_str());

  return vector;
}

std::vector<vtkSmartPointer<Slicer>> fMainWindow::GetComparisonViewers()
{
  std::vector<vtkSmartPointer<Slicer>>comparisonViewers;
  if (mSlicerManagers.size() == 2)
  {
	  comparisonViewers.push_back(m_ComparisonViewerLeft);
	  comparisonViewers.push_back(m_ComparisonViewerCenter);
  }
  else if (mSlicerManagers.size() == 3)
  {
	  comparisonViewers.push_back(m_ComparisonViewerLeft);
	  comparisonViewers.push_back(m_ComparisonViewerCenter);
	  comparisonViewers.push_back(m_ComparisonViewerRight);
  }
  
  return comparisonViewers;
}

void fMainWindow::Registration(std::string fixedFileName, std::vector<std::string> inputFileNames,
  std::vector<std::string> outputFileNames, std::vector<std::string> matrixFileNames,
  std::string metrics, bool rigidMode, bool affineMode, bool deformMode,
  std::string radii, std::string iterations)
{
  std::string configPathName;
  std::string configFileName;
  std::string extn = ".txt";

  std::vector<std::string> affineMatrix;
  std::vector<std::string> outputImage;

  updateProgress(5, "Starting Registration");

  //auto TargetImage = cbica::ReadImage< ImageTypeFloat3D >(fixedFileName);

  if (outputFileNames.size() != inputFileNames.size() || outputFileNames.size() != matrixFileNames.size() || matrixFileNames.size() != inputFileNames.size())
  {
    ShowErrorMessage("Number of input, matrix and output file names do not match");
    return;
  }

  configPathName = itksys::SystemTools::GetFilenamePath(matrixFileNames[0]).c_str();
  configFileName = configPathName + "/" + itksys::SystemTools::GetFilenameWithoutExtension(matrixFileNames[0]).c_str() + extn;

  for (unsigned int i = 0; i < inputFileNames.size(); i++)
  {
    if (!cbica::isFile(inputFileNames[i]))
    {
      ShowErrorMessage("Input file '" + std::to_string(i) + "' is undefined; please check");
      return;
    }
    updateProgress(static_cast<int>(100 / ((i + 1) * inputFileNames.size())), "processing Registration");

    QStringList args;

    args << "-i" << inputFileNames[i].c_str();
    args << "-o" << outputFileNames[i].c_str();
    args << "-rIA" << matrixFileNames[i].c_str();
    args << "-rFI" << fixedFileName.c_str();
    args << "-rNI" << iterations.c_str();

    if (metrics == "NCC")
      args << ("-rME NCC-" + radii).c_str();
    else
      args << "-rME " << metrics.c_str();

    args << "-reg";
    if (rigidMode)
    {
      args << "Rigid";
    }
    else if (affineMode)
    {
      args << "Affine";
    }
    else
    {
      args << "Deformable";
    }
    std::string fullCommandToRun = getApplicationPath("Preprocessing");

    if (startExternalProcess(fullCommandToRun.c_str(), args) != 0)
    {
      ShowErrorMessage("Couldn't register with the default parameters; please use command line functionality");
      return;
    }
    else
    {
      affineMatrix.push_back(matrixFileNames[i] + ".mat");
    }

    if (matrixFileNames[i].find("remove") != std::string::npos)
    {
      if (cbica::isFile(matrixFileNames[i]))
      {
        if (std::remove(matrixFileNames[i].c_str()) == 0)
        {
          updateProgress(80, "Cleaning temporary files");
        }
      }

      updateProgress(static_cast<int>(100 / ((i + 1) * inputFileNames.size())), "Writing File");
    }

    updateProgress(100, "Registration Complete.");

    time_t t = std::time(0);
    long int now = static_cast<long int> (t);

    std::ofstream file;
    file.open(configFileName.c_str());

    std::string mode;

    if (affineMode)
      mode = "Affine";
    else if (rigidMode)
      mode = "Rigid";
    else
      mode = "Deformable";

    if (file.is_open())
    {
      if (metrics != "NCC") {
        file << fixedFileName << ","
          << metrics << ","
          << mode << ","
          << iterations << ","
          << now << "\n";
      }
      else {
        file << fixedFileName << ","
          << metrics << ","
          << radii << ","
          << mode << ","
          << iterations << ","
          << now << "\n";
      }
    }
    file.close();
  }
  //// This happens because the qconcurrent doesn't allow more than 5 function parameters, without std::bind + not sure what else
  //std::vector<std::string> compVector = {
  //  fixedFileName,
  //  ((registrationMode) ? "true" : "false"),
  //  metrics,
  //  ((affineMode) ? "true" : "false"),
  //  radii,
  //  iterations
  //};

  //QtConcurrent::run(this, &fMainWindow::RegistrationWorker,
  //  compVector,
  //  inputFileNames,
  //  outputFileNames,
  //  matrixFileNames
  //);
  /*QFuture<void> r = QtConcurrent::run(std::bind(
    this, &fMainWindow::RegistrationWorker,
    fixedFileName, inputFileNames, outputFileNames,
    matrixFileNames, registrationMode, metrics, affineMode, radii, iterations
  ));*/
}

void fMainWindow::RegistrationWorker(std::vector<std::string> compVector, std::vector<std::string> inputFileNames,
  std::vector<std::string> outputFileNames, std::vector<std::string> matrixFileNames)
{
  // "Unpacking" the variables
  std::string fixedFileName = compVector[0];
  bool registrationMode = (compVector[1] == "true");
  std::string metrics = compVector[2];
  bool affineMode = (compVector[3] == "true");
  std::string radii = compVector[4];
  std::string iterations = compVector[5];

  std::string configPathName;
  std::string configFileName;
  std::string extn = ".txt";

  std::vector<std::string> affineMatrix;
  std::vector<std::string> outputImage;

  updateProgress(5, "Starting Registration");

  //auto TargetImage = cbica::ReadImage< ImageTypeFloat3D >(fixedFileName);

  if (outputFileNames.size() != inputFileNames.size() || outputFileNames.size() != matrixFileNames.size() || matrixFileNames.size() != inputFileNames.size())
  {
    ShowErrorMessage("Number of input, matrix and output file names do not match");
    return;
  }

  configPathName = itksys::SystemTools::GetFilenamePath(matrixFileNames[0]).c_str();
  configFileName = configPathName + "/" + itksys::SystemTools::GetFilenameWithoutExtension(matrixFileNames[0]).c_str() + extn;

  for (unsigned int i = 0; i < inputFileNames.size(); i++)
  {
    if (!cbica::isFile(inputFileNames[i]))
    {
      ShowErrorMessage("Input file '" + std::to_string(i) + "' is undefined; please check");
      return;
    }
    updateProgress(static_cast<int>(100 / ((i + 1) * inputFileNames.size())), "processing Registration");

    std::string fixedFileCommand = "-f " + fixedFileName;
    std::string movingFileCommand = " -i " + inputFileNames[i];
    std::string affineMatrixCommand = " -t " + matrixFileNames[i];
    std::string outputCommand = " -o " + outputFileNames[i];
    std::string metricsCommand = " -m " + metrics;
    std::string iterationsCommand = " -n " + iterations;
    QStringList args;
    args << "-reg" << "-trf" << "-a" << "-f" << fixedFileName.c_str()
      << "-i" << inputFileNames[i].c_str() << "-t" << matrixFileNames[i].c_str() << "-o" << outputFileNames[i].c_str()
      << "-m" << metrics.c_str() << "-n" << iterations.c_str();
        
    if (metrics == "NCC")
      args << "-ri" << radii.c_str();
    if (affineMode)
    {
      args << "-a";
    }
    else
    {
      args << "-r";
    }
    std::string fullCommandToRun = getApplicationPath("GreedyRegistration");

    if (startExternalProcess(fullCommandToRun.c_str(), args) != 0)
    {
      ShowErrorMessage("Couldn't register with the default parameters; please use command line functionality");
      return;
    }
    else
    {
      affineMatrix.push_back(matrixFileNames[i] + ".mat");
    }

    if (matrixFileNames[i].find("remove") != std::string::npos)
    {
      if (cbica::isFile(matrixFileNames[i]))
      {
        if (std::remove(matrixFileNames[i].c_str()) == 0)
        {
          updateProgress(80, "Cleaning temporary files");
        }
      }

      updateProgress(static_cast<int>(100 / ((i + 1) * inputFileNames.size())), "Writing File");
    }

    updateProgress(100, "Registration Complete.");

    time_t t = std::time(0);
    long int now = static_cast<long int> (t);

    std::ofstream file;
    file.open(configFileName.c_str());

    std::string mode;

    if (affineMode == true)
      mode = "Affine";
    else
      mode = "Rigid";

    if (file.is_open())
    {
      if (metrics != "NCC") {
        file << fixedFileName << ","
          << metrics << ","
          << mode << ","
          << iterations << ","
          << now << "\n";
      }
      else {
        file << fixedFileName << ","
          << metrics << ","
          << radii << ","
          << mode << ","
          << iterations << ","
          << now << "\n";
      }
    }
    file.close();
  }

  //std::terminate();
}


void fMainWindow::UpdateAction(std::vector<PointVal> points)
{
  mActionPoints.push_back(points);
}

void fMainWindow::FillLabel(int label)
{
  //auto orientation = mSlicerManagers[0]->mSlicers[0]->GetOrientation();
}

void fMainWindow::UndoFunctionality()
{
  if (mActionPoints.empty())
    return;
  std::vector<PointVal>  OneStrkePoints = mActionPoints.back();
  //Its important to do the undo in reverse order of what happened
  for (std::vector<PointVal>::iterator it = OneStrkePoints.end(); it != OneStrkePoints.begin();)
  {
    --it;
    PointVal pt = *it;
    float* pData = (float*)this->mSlicerManagers[0]->GetSlicer(0)->mMask->GetScalarPointer(pt.x, pt.y, pt.z);
    *pData = pt.value;
  }
  mActionPoints.pop_back();

  this->mSlicerManagers[0]->GetSlicer(0)->mMask->Modified();
  this->mSlicerManagers[0]->Render();
}


void fMainWindow::SetOpacity()
{
  if (this->mSlicerManagers[0]->GetSlicer(0)->GetMaskOpacity() == 0)
    ChangeMaskOpacity(drawingPanel->getCurrentOpacity());
  else
    ChangeMaskOpacity(0);
}

void fMainWindow::closeEvent(QCloseEvent* event)
{
  if (m_NumberOfUnfinishedExternalProcesses > 0)
  {
    ShowErrorMessage("Please close all external applications before exiting.");
    event->ignore();
    return;
  }

  if (m_IsGeodesicTrainingRunning)
  {
    ShowErrorMessage("Please wait for GeodesicTraining execution to finish.");
    event->ignore();
    return;
  }

  if (!cbica::fileExists(closeConfirmation))
  {
    auto msgBox = new QMessageBox();
    msgBox->setWindowTitle("Close Confirmation!");
    msgBox->setText("Are you certain you would like to exit?");
    msgBox->addButton(QMessageBox::Yes);
    msgBox->addButton(QMessageBox::No);
    msgBox->setDefaultButton(QMessageBox::No);

    QCheckBox closeConfirmationBox("Never ask again");
    closeConfirmationBox.blockSignals(true);
    msgBox->addButton(&closeConfirmationBox, QMessageBox::ResetRole);
    if (msgBox->exec() == QMessageBox::Yes)
    {
      if (closeConfirmationBox.checkState() == Qt::Checked)
      {
        std::ofstream file;
        file.open(closeConfirmation.c_str());
        file << "User doesn't want close confirmation.\n";
        file.close();
      }

      //! close the help dialog forcefully as we are about to exit the application
      bool closed = mHelpDlg->close();

      event->accept();
    }
    else
    {
      event->ignore();
    }
  }
  else
  {
    //! close the help dialog forcefully as we are about to exit the application
    bool closed = mHelpDlg->close();

    event->accept();
  }
}

void fMainWindow::updateProgress(int progress, std::string message, int max)
{
#ifdef USE_PROCESSDIALOG
  m_progressBar->setMaximum(max);
  m_progressBar->setValue(progress);
  m_messageLabel->setText(QString::fromStdString(message));
  QTimer::singleShot(10000.0, m_messageLabel, SLOT(clear()));
  qApp->processEvents();

#endif
}
std::vector<std::map<CAPTK::ImageModalityType, std::string>> fMainWindow::LoadQualifiedSubjectsFromGivenDirectoryForSurvival(const std::string directoryname)
{
  std::map<CAPTK::ImageModalityType, std::string> OneQualifiedSubject;
  std::vector<std::map<CAPTK::ImageModalityType, std::string>> QualifiedSubjects;
  std::vector<std::string> subjectNames = cbica::subdirectoriesInDirectory(directoryname);
  std::sort(subjectNames.begin(), subjectNames.end());

  for (unsigned int sid = 0; sid < subjectNames.size(); sid++)
  {
    std::string subjectPath = directoryname + "/" + subjectNames[sid];

    std::string t1ceFilePath = "";
    std::string t1FilePath = "";
    std::string t2FilePath = "";
    std::string t2FlairFilePath = "";
    std::string axFilePath = "";
    std::string faFilePath = "";
    std::string radFilePath = "";
    std::string trFilePath = "";
    std::string rcbvFilePath = "";
    std::string psrFilePath = "";
    std::string phFilePath = "";
    std::string labelPath = "";
    std::string atlasPath = "";
    std::string parametersPath = "";
    std::string featureFilePath = "";

    std::vector<std::string> files;

    if (cbica::directoryExists(subjectPath + "/SEGMENTATION"))
    {
      files = cbica::filesInDirectory(subjectPath + "/SEGMENTATION", false);
      if (files.size() == 1)
      {
        labelPath = subjectPath + "/SEGMENTATION" + "/" + files[0];
      }
      else
      {
        for (unsigned int i = 0; i < files.size(); i++)
        {
          std::string filePath = subjectPath + "/SEGMENTATION" + "/" + files[i], filePath_lower;
          std::string extension = cbica::getFilenameExtension(filePath, false);
          filePath_lower = filePath;
          std::transform(filePath_lower.begin(), filePath_lower.end(), filePath_lower.begin(), ::tolower);
          if ((filePath_lower.find("atlas") != std::string::npos || filePath_lower.find("jakob_label") != std::string::npos)
            && isExtensionSupported(extension))
            atlasPath = subjectPath + "/SEGMENTATION" + "/" + files[i];
          else if ((filePath_lower.find("segmentation") != std::string::npos)
            && isExtensionSupported(extension))
            labelPath = subjectPath + "/SEGMENTATION" + "/" + files[i];
          else if ((filePath_lower.find("parameter") != std::string::npos)
            && (extension == PARAM_EXT))
            parametersPath = subjectPath + "/SEGMENTATION" + "/" + files[i];
        }
      }
    }

    if (cbica::directoryExists(subjectPath + "/CONVENTIONAL"))
    {
      files = cbica::filesInDirectory(subjectPath + "/CONVENTIONAL", false);
      for (unsigned int i = 0; i < files.size(); i++)
      {
        std::string filePath = subjectPath + "/CONVENTIONAL" + "/" + files[i];
        std::string extension = cbica::getFilenameExtension(filePath, false);

        if ((guessImageType(files[i]) == CAPTK::ImageModalityType::IMAGE_TYPE_T1CE) && isExtensionSupported(extension))
          t1ceFilePath = subjectPath + "/CONVENTIONAL" + "/" + files[i];
        else if ((guessImageType(files[i]) == CAPTK::ImageModalityType::IMAGE_TYPE_T1) && isExtensionSupported(extension))
          t1FilePath = subjectPath + "/CONVENTIONAL" + "/" + files[i];
        else if ((guessImageType(files[i]) == CAPTK::ImageModalityType::IMAGE_TYPE_T2) && isExtensionSupported(extension))
          t2FilePath = subjectPath + "/CONVENTIONAL" + "/" + files[i];
        else if ((guessImageType(files[i]) == CAPTK::ImageModalityType::IMAGE_TYPE_T2FLAIR) && isExtensionSupported(extension))
          t2FlairFilePath = subjectPath + "/CONVENTIONAL" + "/" + files[i];
      }
    }

    if (cbica::directoryExists(subjectPath + "/PERFUSION"))
    {
      files = cbica::filesInDirectory(subjectPath + "/PERFUSION", false);
      for (unsigned int i = 0; i < files.size(); i++)
      {
        std::string filePath = subjectPath + "/PERFUSION" + "/" + files[i], filePath_lower;
        std::string extension = cbica::getFilenameExtension(filePath, false);
        filePath_lower = filePath;
        std::transform(filePath_lower.begin(), filePath_lower.end(), filePath_lower.begin(), ::tolower);
        if ((guessImageType(files[i]) == CAPTK::ImageModalityType::IMAGE_TYPE_RCBV)
          && isExtensionSupported(extension))
          rcbvFilePath = subjectPath + "/PERFUSION" + "/" + files[i];
        else if ((guessImageType(files[i]) == CAPTK::ImageModalityType::IMAGE_TYPE_PSR)
          && isExtensionSupported(extension))
          psrFilePath = subjectPath + "/PERFUSION" + "/" + files[i];
        else if ((guessImageType(files[i]) == CAPTK::ImageModalityType::IMAGE_TYPE_PH)
          && isExtensionSupported(extension))
          phFilePath = subjectPath + "/PERFUSION" + "/" + files[i];
      }
    }

    if (cbica::directoryExists(subjectPath + "/DTI"))
    {
      files = cbica::filesInDirectory(subjectPath + "/DTI", false);
      for (unsigned int i = 0; i < files.size(); i++)
      {
        std::string filePath = subjectPath + "/DTI/" + files[i];
        std::string extension = cbica::getFilenameExtension(filePath, false);

        if ((guessImageType(files[i]) == CAPTK::ImageModalityType::IMAGE_TYPE_AX) && isExtensionSupported(extension))
          axFilePath = subjectPath + "/DTI/" + files[i];
        else if ((guessImageType(files[i]) == CAPTK::ImageModalityType::IMAGE_TYPE_FA) && isExtensionSupported(extension))
          faFilePath = subjectPath + "/DTI/" + files[i];
        else if ((guessImageType(files[i]) == CAPTK::ImageModalityType::IMAGE_TYPE_RAD) && isExtensionSupported(extension))
          radFilePath = subjectPath + "/DTI/" + files[i];
        else if ((guessImageType(files[i]) == CAPTK::ImageModalityType::IMAGE_TYPE_TR) && isExtensionSupported(extension))
          trFilePath = subjectPath + "/DTI/" + files[i];
      }
    }
    if (cbica::fileExists(subjectPath + "/features.csv"))
      featureFilePath = subjectPath + "/features.csv";

    if (labelPath.empty() || t1FilePath.empty() || t2FilePath.empty() || t1ceFilePath.empty() || t2FlairFilePath.empty() || rcbvFilePath.empty() || axFilePath.empty() || faFilePath
      == "" || radFilePath.empty() || trFilePath.empty() || psrFilePath.empty() || phFilePath.empty() || featureFilePath.empty())
      continue;

    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_T1] = t1FilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_T2] = t2FilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_T1CE] = t1ceFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_T2FLAIR] = t2FlairFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_AX] = axFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_FA] = faFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_RAD] = radFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_TR] = trFilePath;

    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_PSR] = psrFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_PH] = phFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_RCBV] = rcbvFilePath;

    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_SEG] = labelPath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_ATLAS] = atlasPath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_PARAMS] = parametersPath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_FEATURES] = featureFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_SUDOID] = subjectNames[sid];

    QualifiedSubjects.push_back(OneQualifiedSubject);
  }
  return QualifiedSubjects;
}

std::vector<std::map<CAPTK::ImageModalityType, std::string>> fMainWindow::LoadQualifiedSubjectsFromGivenDirectoryForPCA(const std::string directoryname)
{
  std::map<CAPTK::ImageModalityType, std::string> OneQualifiedSubject;
  std::vector<std::map<CAPTK::ImageModalityType, std::string>> QualifiedSubjects;
  std::vector<std::string> subjectNames = cbica::subdirectoriesInDirectory(directoryname);
  std::sort(subjectNames.begin(), subjectNames.end());

  for (unsigned int sid = 0; sid < subjectNames.size(); sid++)
  {
    std::string subjectPath = directoryname + "/" + subjectNames[sid];

    std::string perfFilePath = "";
    std::string labelPath = "";

    std::vector<std::string> files;
    files = cbica::filesInDirectory(subjectPath + "", false);

    for (unsigned int i = 0; i < files.size(); i++)
    {
      std::string filePath = subjectPath + "/" + files[i], filePath_lower;
      std::string extension = cbica::getFilenameExtension(filePath, false);
      filePath_lower = filePath;
      std::transform(filePath_lower.begin(), filePath_lower.end(), filePath_lower.begin(), ::tolower);
      if ((guessImageType(files[i]) == CAPTK::ImageModalityType::IMAGE_TYPE_SEG)
        && isExtensionSupported(extension))
        labelPath = subjectPath + "/" + files[i];
      else if ((guessImageType(files[i]) == CAPTK::ImageModalityType::IMAGE_TYPE_PERFUSION)
        && isExtensionSupported(extension))
        perfFilePath = subjectPath + "/SEGMENTATION" + "/" + files[i];
    }

    if (labelPath.empty() || perfFilePath.empty())
      continue;

    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_PERFUSION] = perfFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_SEG] = labelPath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_SUDOID] = subjectNames[sid];

    QualifiedSubjects.push_back(OneQualifiedSubject);
  }
  return QualifiedSubjects;
}

std::vector<std::map<CAPTK::ImageModalityType, std::string>>  fMainWindow::LoadQualifiedSubjectsFromGivenDirectoryForRecurrence(const CAPTK::MachineLearningApplicationSubtype type, const std::string &directoryname, const bool &useConventionalData, const bool &useDTIData, const bool &usePerfData, const bool &useDistData)
{
  std::map<CAPTK::ImageModalityType, std::string> OneQualifiedSubject;
  std::vector<std::map<CAPTK::ImageModalityType, std::string>> QualifiedSubjects;
  std::vector<std::string> subjectNames = cbica::subdirectoriesInDirectory(directoryname);
  std::sort(subjectNames.begin(), subjectNames.end());

  for (unsigned int sid = 0; sid < subjectNames.size(); sid++)
  {
    std::string subjectPath = directoryname + "/" + subjectNames[sid];

    std::string t1ceFilePath = "";
    std::string t1FilePath = "";
    std::string t2FilePath = "";
    std::string t2FlairFilePath = "";
    std::string axFilePath = "";
    std::string faFilePath = "";
    std::string radFilePath = "";
    std::string trFilePath = "";
    std::string perfFilePath = "";
    std::string labelPath = "";
    std::string nearFilePath = "";
    std::string farFilePath = "";

    std::vector<std::string> files;


    if (cbica::directoryExists(subjectPath + "/DRAWING"))
    {
      files = cbica::filesInDirectory(subjectPath + "/DRAWING", false);
      for (unsigned int i = 0; i < files.size(); i++)
      {
        std::string filePath = subjectPath + "/DRAWING/" + files[i];
        std::string extension = cbica::getFilenameExtension(filePath, false);
        if ((files[i].find("near") != std::string::npos || files[i].find("Infiltrated") != std::string::npos) && (extension == IMG_EXT || extension == NII_EXT || extension == NII_GZ_EXT))
          nearFilePath = subjectPath + "/DRAWING/" + files[i];

        if ((files[i].find("far") != std::string::npos || files[i].find("Pure") != std::string::npos) && (extension == IMG_EXT || extension == NII_EXT || extension == NII_GZ_EXT))
          farFilePath = subjectPath + "/DRAWING/" + files[i];
      }
    }

    if (cbica::directoryExists(subjectPath + "/SEGMENTATION"))
    {
      files = cbica::filesInDirectory(subjectPath + "/SEGMENTATION", false);
      for (unsigned int i = 0; i < files.size(); i++)
      {
        std::string filePath = subjectPath + "/SEGMENTATION/" + files[i];
        std::string extension = cbica::getFilenameExtension(filePath, false);
        if ((guessImageType(files[i]) == CAPTK::ImageModalityType::IMAGE_TYPE_SEG) && isExtensionSupported(extension))
          labelPath = subjectPath + "/SEGMENTATION/" + files[i];
      }
    }
    if (cbica::directoryExists(subjectPath + "/PERFUSION"))
    {
      files = cbica::filesInDirectory(subjectPath + "/PERFUSION", false);
      for (unsigned int i = 0; i < files.size(); i++)
      {
        std::string filePath = subjectPath + "/PERFUSION/" + files[i];
        std::string extension = cbica::getFilenameExtension(filePath, false);
        if ((guessImageType(files[i]) == CAPTK::ImageModalityType::IMAGE_TYPE_PERFUSION) && isExtensionSupported(extension))
          perfFilePath = subjectPath + "/PERFUSION/" + files[i];
      }
    }
    if (useConventionalData && cbica::directoryExists(subjectPath + "/CONVENTIONAL"))
    {
      files = cbica::filesInDirectory(subjectPath + "/CONVENTIONAL", false);
      for (unsigned int i = 0; i < files.size(); i++)
      {
        std::string filePath = subjectPath + "/CONVENTIONAL/" + files[i];
        std::string extension = cbica::getFilenameExtension(filePath, false);

        if ((guessImageType(files[i]) == CAPTK::ImageModalityType::IMAGE_TYPE_T1CE) && isExtensionSupported(extension))
          t1ceFilePath = subjectPath + "/CONVENTIONAL/" + files[i];
        else if ((guessImageType(files[i]) == CAPTK::ImageModalityType::IMAGE_TYPE_T1) && isExtensionSupported(extension))
          t1FilePath = subjectPath + "/CONVENTIONAL/" + files[i];
        else if ((guessImageType(files[i]) == CAPTK::ImageModalityType::IMAGE_TYPE_T2) && isExtensionSupported(extension))
          t2FilePath = subjectPath + "/CONVENTIONAL/" + files[i];
        else if ((guessImageType(files[i]) == CAPTK::ImageModalityType::IMAGE_TYPE_T2FLAIR) && isExtensionSupported(extension))
          t2FlairFilePath = subjectPath + "/CONVENTIONAL/" + files[i];
      }
    }


    if (useDTIData && cbica::directoryExists(subjectPath + "/DTI"))
    {
      files = cbica::filesInDirectory(subjectPath + "/DTI", false);
      for (unsigned int i = 0; i < files.size(); i++)
      {
        std::string filePath = subjectPath + "/DTI/" + files[i];
        std::string extension = cbica::getFilenameExtension(filePath, false);

        if ((guessImageType(files[i]) == CAPTK::ImageModalityType::IMAGE_TYPE_AX) && isExtensionSupported(extension))
          axFilePath = subjectPath + "/DTI/" + files[i];
        else if ((guessImageType(files[i]) == CAPTK::ImageModalityType::IMAGE_TYPE_FA) && isExtensionSupported(extension))
          faFilePath = subjectPath + "/DTI/" + files[i];
        else if ((guessImageType(files[i]) == CAPTK::ImageModalityType::IMAGE_TYPE_RAD) && isExtensionSupported(extension))
          radFilePath = subjectPath + "/DTI/" + files[i];
        else if ((guessImageType(files[i]) == CAPTK::ImageModalityType::IMAGE_TYPE_TR) && isExtensionSupported(extension))
          trFilePath = subjectPath + "/DTI/" + files[i];
      }
    }


    if ((useConventionalData && t1FilePath.empty()) || (useConventionalData && t2FilePath.empty()) || (useConventionalData && t1ceFilePath.empty()) || (useConventionalData && t2FlairFilePath.empty()) ||
      (usePerfData && perfFilePath.empty()) || (useDTIData && axFilePath.empty()) || (useDTIData && faFilePath.empty()) || (useDTIData && radFilePath.empty()) || (useDTIData && trFilePath.empty()))
      continue;

    if ((nearFilePath.empty() || farFilePath.empty()) && type == CAPTK::MachineLearningApplicationSubtype::TRAINING)
      continue;

    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_T1] = t1FilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_T2] = t2FilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_T1CE] = t1ceFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_T2FLAIR] = t2FlairFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_AX] = axFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_FA] = faFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_RAD] = radFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_TR] = trFilePath;

    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_SEG] = labelPath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_SUDOID] = subjectNames[sid];
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_PERFUSION] = perfFilePath;

    if (type == CAPTK::MachineLearningApplicationSubtype::TRAINING)
    {
      OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_NEAR] = nearFilePath;
      OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_FAR] = farFilePath;
    }
    QualifiedSubjects.push_back(OneQualifiedSubject);

  }
  return QualifiedSubjects;
}


std::vector<std::map<CAPTK::ImageModalityType, std::string>>  fMainWindow::LoadQualifiedSubjectsFromGivenDirectoryForPseudoProgression(const CAPTK::MachineLearningApplicationSubtype type, const std::string &directoryname, const bool &useConventionalData, const bool &useDTIData, const bool &usePerfData, const bool &useDistData)
{
  std::map<CAPTK::ImageModalityType, std::string> OneQualifiedSubject;
  std::vector<std::map<CAPTK::ImageModalityType, std::string>> QualifiedSubjects;
  std::vector<std::string> subjectNames = cbica::subdirectoriesInDirectory(directoryname);

  //subjectNames.clear();
  //subjectNames.push_back("AAMA");
  //subjectNames.push_back("AAMG");
  //subjectNames.push_back("AAMJ");
  //subjectNames.push_back("AAMP");
  //subjectNames.push_back("AAMQ");
  //subjectNames.push_back("ABEM");

  std::sort(subjectNames.begin(), subjectNames.end());

  for (unsigned int sid = 0; sid < subjectNames.size(); sid++)
  {
    std::string subjectPath = directoryname + "/" + subjectNames[sid];

    std::string t1ceFilePath = "";
    std::string t1FilePath = "";
    std::string t2FilePath = "";
    std::string t2FlairFilePath = "";

    //std::string t1t1ceFilePath    = "";
    //std::string t2t2FlairFilePath = "";

    std::string axFilePath = "";
    std::string faFilePath = "";
    std::string radFilePath = "";
    std::string trFilePath = "";
    std::string phFilePath = "";
    std::string psrFilePath = "";
    std::string rcbvFilePath = "";
    std::string perfFilePath = "";
    std::string featuresFilePath = "";

    std::string labelPath = "";
    std::string atlasPath = "";

    std::vector<std::string> files;


    if (cbica::directoryExists(subjectPath + "/SEGMENTATION"))
    {
      files = cbica::filesInDirectory(subjectPath + "/SEGMENTATION", false);
      for (unsigned int i = 0; i < files.size(); i++)
      {
        std::string filePath = subjectPath + "/SEGMENTATION/" + files[i];
        std::string extension = cbica::getFilenameExtension(filePath, false);
        if ((guessImageType(files[i]) == CAPTK::ImageModalityType::IMAGE_TYPE_SEG) && isExtensionSupported(extension))
          labelPath = subjectPath + "/SEGMENTATION/" + files[i];
        else if ((files[i].find("atlas") != std::string::npos) && isExtensionSupported(extension))
          atlasPath = subjectPath + "/SEGMENTATION/" + files[i];
      }
    }
    if (cbica::directoryExists(subjectPath + "/PERFUSION"))
    {
      files = cbica::filesInDirectory(subjectPath + "/PERFUSION", false);
      for (unsigned int i = 0; i < files.size(); i++)
      {
        std::string filePath = subjectPath + "/PERFUSION" + "/" + files[i], filePath_lower;
        std::string extension = cbica::getFilenameExtension(filePath, false);
        filePath_lower = filePath;
        std::transform(filePath_lower.begin(), filePath_lower.end(), filePath_lower.begin(), ::tolower);
        if ((guessImageType(files[i]) == CAPTK::ImageModalityType::IMAGE_TYPE_RCBV)
          && isExtensionSupported(extension))
          rcbvFilePath = subjectPath + "/PERFUSION" + "/" + files[i];
        else if ((guessImageType(files[i]) == CAPTK::ImageModalityType::IMAGE_TYPE_PSR)
          && isExtensionSupported(extension))
          psrFilePath = subjectPath + "/PERFUSION" + "/" + files[i];
        else if ((guessImageType(files[i]) == CAPTK::ImageModalityType::IMAGE_TYPE_PH)
          && isExtensionSupported(extension))
          phFilePath = subjectPath + "/PERFUSION" + "/" + files[i];
        else if ((guessImageType(files[i]) == CAPTK::ImageModalityType::IMAGE_TYPE_PERFUSION)
          && isExtensionSupported(extension))
          perfFilePath = subjectPath + "/PERFUSION" + "/" + files[i];
      }
    }
    if (useConventionalData && cbica::directoryExists(subjectPath + "/CONVENTIONAL"))
    {
      files = cbica::filesInDirectory(subjectPath + "/CONVENTIONAL", false);
      for (unsigned int i = 0; i < files.size(); i++)
      {
        std::string filePath = subjectPath + "/CONVENTIONAL/" + files[i];
        std::string extension = cbica::getFilenameExtension(filePath, false);

        if ((guessImageType(files[i]) == CAPTK::ImageModalityType::IMAGE_TYPE_T1CE) && isExtensionSupported(extension))
          t1ceFilePath = subjectPath + "/CONVENTIONAL/" + files[i];
        else if ((guessImageType(files[i]) == CAPTK::ImageModalityType::IMAGE_TYPE_T1) && isExtensionSupported(extension))
          t1FilePath = subjectPath + "/CONVENTIONAL/" + files[i];
        else if ((guessImageType(files[i]) == CAPTK::ImageModalityType::IMAGE_TYPE_T2) && isExtensionSupported(extension))
          t2FilePath = subjectPath + "/CONVENTIONAL/" + files[i];
        else if ((guessImageType(files[i]) == CAPTK::ImageModalityType::IMAGE_TYPE_T2FLAIR) && isExtensionSupported(extension))
          t2FlairFilePath = subjectPath + "/CONVENTIONAL/" + files[i];
      }
    }


    if (useDTIData && cbica::directoryExists(subjectPath + "/DTI"))
    {
      files = cbica::filesInDirectory(subjectPath + "/DTI", false);
      for (unsigned int i = 0; i < files.size(); i++)
      {
        std::string filePath = subjectPath + "/DTI/" + files[i];
        std::string extension = cbica::getFilenameExtension(filePath, false);

        if ((guessImageType(files[i]) == CAPTK::ImageModalityType::IMAGE_TYPE_AX) && isExtensionSupported(extension))
          axFilePath = subjectPath + "/DTI/" + files[i];
        else if ((guessImageType(files[i]) == CAPTK::ImageModalityType::IMAGE_TYPE_FA) && isExtensionSupported(extension))
          faFilePath = subjectPath + "/DTI/" + files[i];
        else if ((guessImageType(files[i]) == CAPTK::ImageModalityType::IMAGE_TYPE_RAD) && isExtensionSupported(extension))
          radFilePath = subjectPath + "/DTI/" + files[i];
        else if ((guessImageType(files[i]) == CAPTK::ImageModalityType::IMAGE_TYPE_TR) && isExtensionSupported(extension))
          trFilePath = subjectPath + "/DTI/" + files[i];
      }
    }
    if (cbica::fileExists(subjectPath + "/features.csv"))
      featuresFilePath = subjectPath + "/features.csv";

    if (labelPath.empty() || t1FilePath.empty() || t2FilePath.empty() || t1ceFilePath.empty() || t2FlairFilePath.empty() || rcbvFilePath.empty() || axFilePath.empty() || faFilePath.empty() || radFilePath.empty() || trFilePath.empty() || psrFilePath.empty() || phFilePath.empty())
      continue;

    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_T1] = t1FilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_T2] = t2FilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_T1CE] = t1ceFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_T2FLAIR] = t2FlairFilePath;
    //OneQualifiedSubject[IMAGE_TYPE_T1T1CE]    = t1t1ceFilePath; 
    //OneQualifiedSubject[IMAGE_TYPE_T2FL]      = t2t2FlairFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_AX] = axFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_FA] = faFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_RAD] = radFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_TR] = trFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_PH] = phFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_PSR] = psrFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_RCBV] = rcbvFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_PERFUSION] = perfFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_FEATURES] = featuresFilePath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_SEG] = labelPath;
    OneQualifiedSubject[CAPTK::ImageModalityType::IMAGE_TYPE_SUDOID] = subjectNames[sid];
    QualifiedSubjects.push_back(OneQualifiedSubject);

    std::cout << subjectNames[sid] << std::endl;
  }
  return QualifiedSubjects;
}

bool fMainWindow::isMaskDefined()
{
  if (!mSlicerManagers.empty())
  {
    auto minMaxComputer = itk::MinimumMaximumImageCalculator<ImageTypeFloat3D>::New();
    minMaxComputer->SetImage(convertVtkToItk<float, 3>(mSlicerManagers[0]->mMask));
    minMaxComputer->Compute();
    if (minMaxComputer->GetMaximum() > 0)
    {
      return true;
    }
  }

  return false;
}

void fMainWindow::SetComparisonMode(bool mode)
{
  this->m_ComparisonMode = mode;
}

bool fMainWindow::GetComparisonMode()
{
  return this->m_ComparisonMode;
}

void fMainWindow::help_contextual(const std::string startPage)
{
  mHelpDlg->setNewStartPage(startPage);
  mHelpDlg->show();
}

std::vector< fMainWindow::ActionAndName >fMainWindow::populateStringListInMenu(const std::string &inputList, QMainWindow* inputFMainWindow, QMenu* menuToPopulate, std::string menuAppSubGroup, bool ExcludeGeodesic)
{
  if (!inputList.empty())
  {
    std::string inputList_wrap = inputList;
    if (inputList_wrap[0] == ' ')
    {
      inputList_wrap.erase(0, 1);
    }
    std::vector< std::string > vectorOfInputs = cbica::stringSplit(inputList_wrap, " ");
    return populateStringListInMenu(vectorOfInputs, inputFMainWindow, menuToPopulate, menuAppSubGroup, ExcludeGeodesic);
  }
  else
  {
    return std::vector< fMainWindow::ActionAndName >{};
  }
}

std::vector< fMainWindow::ActionAndName >fMainWindow::populateStringListInMenu(const std::vector< std::string > &vectorOfInputs, QMainWindow* inputFMainWindow, QMenu* menuToPopulate, std::string menuAppSubGroup, bool ExcludeGeodesic)
{
  std::vector< fMainWindow::ActionAndName > returnVector;
  if (ExcludeGeodesic)
  {
    returnVector.resize(vectorOfInputs.size());
  }
  else
  {
    returnVector.resize(vectorOfInputs.size() + 1);
  }
  size_t returnVecCounter = 0;

  if (!menuAppSubGroup.empty())
  {
    returnVector[returnVecCounter].action = new QAction(inputFMainWindow);
    returnVector[returnVecCounter].action->setObjectName(QString::fromUtf8(std::string("action" + menuAppSubGroup).c_str()));
    returnVector[returnVecCounter].action->setIconText(QString(menuAppSubGroup.c_str()));
    returnVector[returnVecCounter].action->setText(QString(menuAppSubGroup.c_str()));
    returnVector[returnVecCounter].action->setEnabled(false);
    returnVector[returnVecCounter].name = menuAppSubGroup;
    returnVector[returnVecCounter].action->setEnabled(false);
    menuToPopulate->addAction(returnVector[returnVecCounter].action);
    returnVecCounter++;
  }

  for (size_t i = 0; i < vectorOfInputs.size(); i++)
  {
    if ((vectorOfInputs[i] != "FeatureExtraction"))
    {
      returnVector[returnVecCounter].action = new QAction(inputFMainWindow);
      returnVector[returnVecCounter].action->setObjectName(QString::fromUtf8(std::string("action" + vectorOfInputs[i]).c_str()));
      returnVector[returnVecCounter].action->setIconText(QString(vectorOfInputs[i].c_str()));
      returnVector[returnVecCounter].name = vectorOfInputs[i];
#ifdef CAPTK_BUILD_CONSOLE_ONLY
      returnVector[returnVecCounter].action->setEnabled(false);
#endif
      menuToPopulate->addAction(returnVector[returnVecCounter].action);
      returnVecCounter++;
    }
  }
  //}

  return returnVector;
}

void fMainWindow::OnPreferencesMenuClicked()
{
	int result = this->preferenceDialog->exec();
}

void fMainWindow::OnSegmentationClicked()
{
  ShowMessage("Inference through the graphical interface is disabled for this release, please try from the command line.");
  return;

  auto outputDir = segmentationPanel->getOutputPath();
  auto inputDir = segmentationPanel->getInputDirectoryPath();
  auto allSegmentationsOptions = segmentationPanel->getSelectedSegmentationAlgorithms();
  auto allFusionOptions = segmentationPanel->getSelectedFusionAlgorithms();
  auto currentModel = segmentationPanel->getSelectedModel();

  if (!cbica::isDir(inputDir))
  {
    ShowErrorMessage("Specified input directory is not present.");
    return;
  }

  if (allSegmentationsOptions.empty())
  {
    ShowErrorMessage("No models have been selected for inference. Please select at least 1 and re-try.");
    return;
  }

  if (allSegmentationsOptions.size() > 1)
  {
    if (allFusionOptions.empty())
    {
      ShowErrorMessage("No fusion methods have been selected. If there are more than 1 models selected for inference, at least 1 fusion method is needed.");
      return;
    }
  }

  for (size_t i = 0; i < allSegmentationsOptions.size(); i++)
  {
    if (allSegmentationsOptions[i] == "DeepMedic")
    {
      auto deepMedicExe = getApplicationPath("DeepMedic");

      auto subDirs = cbica::subdirectoriesInDirectory(inputDir);

      std::string subjectsWithMissingModalities, subjectsWithErrors;

      for (size_t s = 0; s < subDirs.size(); s++)
      {
        std::string file_t1gd, file_t1, file_t2, file_flair;
        auto fileToCheck = inputDir + "/" + subDirs[s] + "/brain_t1gd.nii.gz";
        if (cbica::fileExists(fileToCheck))
        {
          file_t1gd = fileToCheck;
        }
        else
        {
          subjectsWithMissingModalities += subDirs[s];
        }

        fileToCheck = inputDir + "/" + subDirs[s] + "/brain_t1.nii.gz";
        if (cbica::fileExists(fileToCheck))
        {
          file_t1 = fileToCheck;
        }
        else
        {
          subjectsWithMissingModalities += subDirs[s];
        }
        fileToCheck = inputDir + "/" + subDirs[s] + "/brain_t2.nii.gz";
        if (cbica::fileExists(fileToCheck))
        {
          file_t2 = fileToCheck;
        }
        else
        {
          subjectsWithMissingModalities += subDirs[s];
        }
        fileToCheck = inputDir + "/" + subDirs[s] + "/brain_flair.nii.gz";
        if (cbica::fileExists(fileToCheck))
        {
          file_flair = fileToCheck;
        }
        else
        {
          subjectsWithMissingModalities += subDirs[s] + ",";
        }

        if (subjectsWithMissingModalities.empty())
        {
          auto brainMaskFile = inputDir + "/" + subDirs[s] + "/deepmedic_seg.nii.gz";

          auto fullCommand = deepMedicExe + " -md " + getCaPTkDataDir() + "/fets/deepMedic/saved_models/brainTumorSegmentation/ " +
          "-i " + file_t1 + "," +
          file_t1gd + "," +
          file_t2 + "," +
          file_flair + " -o " +
          brainMaskFile;
          
          if (std::system(fullCommand.c_str()) != 0)
          {
            subjectsWithErrors += subDirs[s] + ",";
          }
        }
      }

      if (!subjectsWithMissingModalities.empty())
      {
        subjectsWithMissingModalities.pop_back();
        ShowErrorMessage("The following subjects had at least one missing modality: \n" + subjectsWithMissingModalities);
      }
      if (!subjectsWithErrors.empty())
      {
        subjectsWithErrors.pop_back();
        ShowErrorMessage("DeepMedic couldn't run the following subjects: \n" + subjectsWithErrors);
      }

    }
  }
  
  QString hardcodedPlanName,
    hardcodedModelWeightPath = (captk_currentApplicationPath + "/OpenFederatedLearning/bin/federations/weights/").c_str(), // start with the common location
    hardcodedPythonPath = (captk_currentApplicationPath + "/OpenFederatedLearning/venv/bin/python").c_str(), // this needs to change for Windows (wonder what happens for macOS?)
    hardcodedModelName = "";

  QStringList args;

  if (currentModel == 0)
  {
    hardcodedPlanName = "pt_3dresunet_ss_brainmagebrats";
    hardcodedModelName = hardcodedModelWeightPath + hardcodedPlanName + "_best.pt"; // taken from https://github.com/FETS-AI/Models/blob/master/skullstripping/3dresunet/pt_3dresunet_ss_brainmagebrats_best.pt
    args << "-nmwf" << hardcodedModelName;
  }
  else
  {
    hardcodedPlanName = "pt_3dresunet_brainmagebrats";
    auto hardcodedModelName = hardcodedPlanName + "_best.pbuf";
    if (!cbica::isFile((hardcodedModelWeightPath + "/" + hardcodedModelName).toStdString()))
    {
      auto hardcodedModelName = hardcodedPlanName + "_init.pbuf";
      if (!cbica::isFile((hardcodedModelWeightPath + "/" + hardcodedModelName).toStdString()))
      {
        ShowErrorMessage("A compatible model weight file was not found. Please contact admin@fets.ai for help.");
        return;
      }
    }
    args << "-mwf" << hardcodedModelName;
  }

  // sanity checks
  //if (!cbica::isFile(hardcodedModelWeightPath.toStdString())) // todo: renable after model weights are using full paths for all
  //{
  //  ShowErrorMessage("The requested inference model was not found (it needs to be in ${FeTS_installDir}/bin/OpenFederatedLearning/bin/federations/weights/${planName}_best.pbuf");
  //  return;
  //}
  if (!cbica::isFile(hardcodedPythonPath.toStdString()))
  {
    ShowErrorMessage("The python virtual environment was not found, please refer to documentation to initialize it.");
    return;
  }

  QString fullCommandToRun = (hardcodedPythonPath.toStdString() + " " +
    captk_currentApplicationPath + "/OpenFederatedLearning/bin/run_inference_from_flplan.py").c_str();

  args << "-p" << hardcodedPlanName + ".yaml"
    //<< "-mwf" << hardcodedModelWeightPath // todo: doing customized solution above - change after model weights are using full paths for all
    << "-d" << inputDir.c_str()
    << "-ld" << m_tempFolderLocation.c_str();

  args << "-md";
  if (segmentationPanel->isGPUenabled())
  {
    args << "cuda";
  }
  else
  {
    args << "cpu";
  }

  QMessageBox *box = new QMessageBox(QMessageBox::Question, "Long running Application",
    "Once inference starts, FeTS UI will not be responsive; press OK to continue...",
    QMessageBox::Ok | QMessageBox::Cancel);
  box->setAttribute(Qt::WA_DeleteOnClose); //makes sure the msgbox is deleted automatically when closed
  box->setWindowModality(Qt::NonModal);
  QCoreApplication::processEvents();
  if (box->exec() == QMessageBox::Ok)
  {
    if (startExternalProcess(fullCommandToRun, args) != 0)
    {
      ShowErrorMessage("Couldn't complete the inference task from GUI; please use command line functionality");
      return;
    }
  }
}

void fMainWindow::OnTrainingClicked()
{
  ShowMessage("Training is disabled for this release.");
  return;

  auto outputDir = trainingPanel->getOutputPath();
  auto inputDir = trainingPanel->getInputDirectoryPath();
  auto allSegmentationsOptions = trainingPanel->getSelectedSegmentationAlgorithms();
  auto colCommonName = trainingPanel->getCollaboratorName();

  if (!cbica::isDir(inputDir))
  {
    ShowErrorMessage("Specified input directory is not present.");
    return;
  }

  //// todo: uncomment this when more algorithms are integrated
  //if (allSegmentationsOptions.first.empty())
  //{
  //  ShowErrorMessage("No models have been selected for inference. Please select at least 1 and re-try.");
  //  return;
  //}

  //if (allSegmentationsOptions.second.empty()) // this is a fresh training cycle
  //{
  //  ShowMessage("Fresh Training Cycle selected with Architecture " + allSegmentationsOptions.first);
  //}
  //else
  //{
  //  ShowMessage("Training Cycle selected with Architecture '" + allSegmentationsOptions.first + 
  //    "' and initial configuration from '" + allSegmentationsOptions.second + "'");
  //  // take the initial model configration from 'allSegmentationsOptions.second'
  //}

  // add ui control for gpu/cpu

  QString hardcodedPlanName = "pt_3dresunet_brainmagebrats", // start with the common location
    hardcodedPythonPath = (captk_currentApplicationPath + "/OpenFederatedLearning/venv/bin/python").c_str(); // this needs to change for Windows (wonder what happens for macOS?)

  // sanity checks
  if (!cbica::isFile(hardcodedPythonPath.toStdString()))
  {
    ShowErrorMessage("The python virtual environment was not found, please refer to documentation to initialize it.");
    return;
  }

  QString fullCommandToRun = (hardcodedPythonPath.toStdString() + " " +
    captk_currentApplicationPath + "/OpenFederatedLearning/bin/run_collaborator_from_flplan.py").c_str();
  QStringList args;
  args << "-p" << hardcodedPlanName + ".yaml"
    << "-col" << colCommonName.c_str() // get this from the UI and needs to be match the name in the certificate
    << "-d" << inputDir.c_str()
    << "-ld" << m_tempFolderLocation.c_str();

  args << "-md";
  if (trainingPanel->isGPUenabled())
  {
    args << "cuda";
  }
  else
  {
    args << "cpu";
  }

  QMessageBox *box = new QMessageBox(QMessageBox::Question, "Long running Application",
    "Once training starts, FeTS UI will not be responsive; press OK to continue...",
    QMessageBox::Ok | QMessageBox::Cancel);
  box->setAttribute(Qt::WA_DeleteOnClose); //makes sure the msgbox is deleted automatically when closed
  box->setWindowModality(Qt::NonModal);
  QCoreApplication::processEvents();
  if (box->exec() == QMessageBox::Ok)
  {
    if (startExternalProcess(fullCommandToRun, args) != 0)
    {
      ShowErrorMessage("Couldn't complete the training task from GUI; please use command line functionality");
      return;
    }
  }
}
