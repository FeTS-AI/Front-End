#include "fTrainingPanel.h"
#include "yaml-cpp/yaml.h"

fTrainingPanel::fTrainingPanel(QWidget * parent) : QWidget(parent)
{
  setupUi(this);

  connect(HelpButton, SIGNAL(clicked()), this, SLOT(helpClicked()));
  connect(m_btnCompute, SIGNAL(clicked()), this, SLOT(onComputeButtonClicked()));
  //connect(m_btnBrowseSaveFile, SIGNAL(clicked()), this, SLOT(browseOutputFileName()));
  connect(m_btnInput, SIGNAL(clicked()), this, SLOT(browseInputDirectory()));

  QSignalMapper *segmentationAlgosMapper = new QSignalMapper();

  // each checkbox is connected to their own dialog implementation and is then tied with the output structure that is to be passed on to fMainWindow later
  for (size_t i = 0; i < m_algoRadioButtonMap.size(); i++)
  {
    connect(m_algoRadioButtonMap[i], SIGNAL(toggledOnActive()), segmentationAlgosMapper, SLOT(map()));
    segmentationAlgosMapper->setMapping(m_algoRadioButtonMap[i], i);
  }
  connect(segmentationAlgosMapper, SIGNAL(mapped(int)), this, SLOT(selectTrainingAlgorithm(int)));
}

void fTrainingPanel::selectTrainingAlgorithm(int selectedAlgorithm)
{
  auto algorithmPathToCheck = (m_fetsDataDir + "/" + m_algoRadioButtonMap[selectedAlgorithm]->text()).toStdString();
  if (cbica::isDir(algorithmPathToCheck))
  {
    auto filesInSubDir = cbica::filesInDirectory(algorithmPathToCheck, false);

    auto pathToArch = (m_fetsDataDir + "/" + m_algoRadioButtonMap[selectedAlgorithm]->text()).toStdString();

    // ensure that the previously trained model file is detected ONLY when a compatible configuration is present
    // this check is done in the next if-loop
    m_trainedModelConfigFile.clear(); 

    // construct pop-up only if there are trained configs present
    if (!filesInSubDir.empty())
    {
      m_dialog = new QDialog();

      auto trainedModelsLayout = new QGridLayout();

      // setup the display headers
      auto checkbox_header = new QLabel(" ");
      auto subjects_header = new QLabel("<b>No. of Subjects</b>");
      auto timestamp_header = new QLabel("<b>Timestamp</b>");
      auto architecture_header = new QLabel("<b>Architecture</b>");
      auto filename_header = new QLabel("<b>Filename</b>");

      trainedModelsLayout->addWidget(checkbox_header, 0, 0, Qt::AlignCenter);
      trainedModelsLayout->addWidget(subjects_header, 0, 1, Qt::AlignCenter);
      trainedModelsLayout->addWidget(timestamp_header, 0, 2, Qt::AlignCenter);
      trainedModelsLayout->addWidget(architecture_header, 0, 3, Qt::AlignCenter);
      trainedModelsLayout->addWidget(filename_header, 0, 4, Qt::AlignCenter);

      m_dialog->setWindowTitle(m_algoRadioButtonMap[selectedAlgorithm]->text() + " Models:");

      QString prevArch;
      for (size_t i = 0; i < filesInSubDir.size(); i++)
      {
        m_algoNestedRadioButtonMap[selectedAlgorithm][i] = new QRadioButton(" ", m_dialog);
        m_algoNestedRadioButtonMap[selectedAlgorithm][i]->setChecked(false);
        m_trainedModelConfigFile[selectedAlgorithm][i] = pathToArch + "/" + filesInSubDir[i]; // store the compatible trained model file for re-training

        // metadata
        auto config = YAML::LoadFile(m_trainedModelConfigFile[selectedAlgorithm][i]);

        auto temp = config["subjects"].as< std::string >();

        auto subjects = new QLabel(config["subjects"].as< std::string >().c_str());
        auto timestamp = new QLabel(config["timestamp"].as< std::string >().c_str());
        auto architecture = new QLabel(config["architecture"].as< std::string >().c_str());
        auto filename = new QLabel(filesInSubDir[i].c_str());

        prevArch = architecture->text();

        // put things in order
        trainedModelsLayout->addWidget(m_algoNestedRadioButtonMap[selectedAlgorithm][i], i + 1, 0, Qt::AlignCenter);
        trainedModelsLayout->addWidget(subjects, i + 1, 1, Qt::AlignCenter);
        trainedModelsLayout->addWidget(timestamp, i + 1, 2, Qt::AlignCenter);
        trainedModelsLayout->addWidget(architecture, i + 1, 3, Qt::AlignCenter);
        trainedModelsLayout->addWidget(filename, i + 1, 4, Qt::AlignCenter);
      }

      // this is for a fresh training
      {
        auto freshModelCounter = filesInSubDir.size();
        m_algoNestedRadioButtonMap[selectedAlgorithm][freshModelCounter] = new QRadioButton(" ", m_dialog);
        m_algoNestedRadioButtonMap[selectedAlgorithm][freshModelCounter]->setChecked(false);
        trainedModelsLayout->addWidget(m_algoNestedRadioButtonMap[selectedAlgorithm][freshModelCounter], freshModelCounter + 1, 0, Qt::AlignCenter);

        auto subjects = new QLabel("N.A.");
        auto timestamp = new QLabel("N.A.");
        auto architecture = new QLabel(prevArch);
        auto filename = new QLabel("New-Model");

        trainedModelsLayout->addWidget(subjects, freshModelCounter + 1, 1, Qt::AlignCenter);
        trainedModelsLayout->addWidget(timestamp, freshModelCounter + 1, 2, Qt::AlignCenter);
        trainedModelsLayout->addWidget(architecture, freshModelCounter + 1, 3, Qt::AlignCenter);
        trainedModelsLayout->addWidget(filename, freshModelCounter + 1, 4, Qt::AlignCenter);
      }

      auto buttonRowCounter = filesInSubDir.size() + 1;

      auto buttonWidth = QLabel().fontMetrics().width("ButtonSize------------------");
      auto buttonHeight = QLabel().fontMetrics().width("----------");

      // initialize the ok and cancel buttons and make them do something meaningful
      auto okButton = new QPushButton("Ok");
      auto cancelButton = new QPushButton("Cancel");

      okButton->setFixedWidth(buttonHeight * 2);
      cancelButton->setFixedWidth(buttonHeight * 2);
      trainedModelsLayout->addWidget(okButton, buttonRowCounter + 1, 3, Qt::AlignCenter);
      trainedModelsLayout->addWidget(cancelButton, buttonRowCounter + 1, 4, Qt::AlignCenter);
      connect(okButton, SIGNAL(clicked()), this, SLOT(onOkButtonClicked()));
      connect(cancelButton, SIGNAL(clicked()), this, SLOT(onCancelButtonClicked()));

      m_dialog->setMaximumSize(buttonWidth * 4, buttonHeight * buttonRowCounter);
      m_dialog->setLayout(trainedModelsLayout);
      m_dialog->setWindowModality(Qt::ApplicationModal);
      m_dialog->exec();
    }
  }
}

void fTrainingPanel::onOkButtonClicked()
{
  // no need to do anything here; all control is happening in getSelectedSegmentationAlgorithms()
  m_dialog->close();
}

void fTrainingPanel::onCancelButtonClicked()
{
  // no untoggle is needed since RadioButtons do that automatically
  m_dialog->close();
}

//void fTrainingPanel::browseOutputFileName()
//{
//  m_outputPath = getExistingDirectory(this, mInputPathName);
//}

void fTrainingPanel::browseInputDirectory()
{
  m_inputDataDirectory = getExistingDirectory(this, mInputPathName);
  m_txtInputDirName->setText(m_inputDataDirectory);
  m_txtSaveFileName->setText(m_inputDataDirectory);
}

void fTrainingPanel::helpClicked()
{
  emit helpClicked_SegmentationUsage("gs_drawing.html");
}

void fTrainingPanel::onComputeButtonClicked()
{
  emit m_btnComputeClicked();
}

std::pair< std::string, std::string > fTrainingPanel::getSelectedSegmentationAlgorithms()
{
  std::string selectedArch;
  for (size_t i = 0; i < m_algoRadioButtonMap.size(); i++) // loop through all the segmentation algorithms
  {
    if (m_algoRadioButtonMap[i]->isChecked())
    {
      selectedArch = m_algoRadioButtonMap[i]->text().toStdString();
      for (size_t j = 0; j < m_algoNestedRadioButtonMap[i].size(); j++)
      {
        if (m_algoNestedRadioButtonMap[i][j]->isChecked())
        {
          return std::make_pair(selectedArch, m_trainedModelConfigFile[i][j]);
        }
      }
    }
  }

  // no prior trained config found for the selected architecture
  return std::make_pair(selectedArch, "");
}
