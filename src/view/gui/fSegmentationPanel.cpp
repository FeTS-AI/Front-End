#include "fSegmentationPanel.h"
#include "yaml-cpp/yaml.h"

#include "qsignalmapper.h"

fSegmentationPanel::fSegmentationPanel(QWidget * parent) : QWidget(parent)
{
  setupUi(this);

  connect(HelpButton, SIGNAL(clicked()), this, SLOT(helpClicked()));
  connect(m_btnCompute, SIGNAL(clicked()), this, SLOT(onComputeButtonClicked()));
  //connect(m_btnBrowseSaveFile, SIGNAL(clicked()), this, SLOT(browseOutputFileName()));
  connect(m_btnInput, SIGNAL(clicked()), this, SLOT(browseInputDirectory()));
  
  QSignalMapper *segmentationAlgosMapper = new QSignalMapper();
  QSignalMapper *segmentationResetAlgosMapper = new QSignalMapper();

  // each checkbox is connected to their own dialog implementation and is then tied with the output structure that is to be passed on to fMainWindow later
  for (auto it = m_algoCheckBoxMap.begin(); it != m_algoCheckBoxMap.end(); ++it)
  {
    connect(it->second, SIGNAL(clickedOnActive()), segmentationAlgosMapper, SLOT(map()));
    connect(it->second, SIGNAL(clickedOnInActive()), segmentationResetAlgosMapper, SLOT(map()));
    segmentationAlgosMapper->setMapping(it->second, it->first);
    segmentationResetAlgosMapper->setMapping(it->second, it->first);
  }
  connect(segmentationAlgosMapper, SIGNAL(mapped(int)), this, SLOT(selectSegmentationAlgorithm(int)));
  connect(segmentationResetAlgosMapper, SIGNAL(mapped(int)), this, SLOT(unSelectSegmentationAlgorithm(int)));
}

//void fSegmentationPanel::browseOutputFileName()
//{
//  m_outputPath = getExistingDirectory(this, mInputPathName);
//}

void fSegmentationPanel::browseInputDirectory()
{
  m_inputDataDirectory = getExistingDirectory(this, mInputPathName);
  m_txtInputDirName->setText(m_inputDataDirectory);
  m_txtSaveFileName->setText(m_inputDataDirectory);
}

void fSegmentationPanel::unSelectSegmentationAlgorithm(int selectedAlgorithm)
{
  m_algoCounterMap[m_currentSegmentationAlgorithm]->setText("0");
}

void fSegmentationPanel::selectSegmentationAlgorithm(int selectedAlgorithm)
{
  m_currentSegmentationAlgorithm = selectedAlgorithm;

  auto pathToArch = (m_fetsDataDir + "/" + m_algoCheckBoxMap[selectedAlgorithm]->text()).toStdString();

  auto filesInSubDir = cbica::filesInDirectory(pathToArch, false);

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

    m_dialog->setWindowTitle(m_algoCheckBoxMap[selectedAlgorithm]->text() + " Models:");

    for (size_t i = 0; i < filesInSubDir.size(); i++)
    {
      m_algoNestedCheckBoxMap[selectedAlgorithm][i] = new QCheckBox(" ", m_dialog);
      m_algoNestedCheckBoxMap[selectedAlgorithm][i]->setChecked(false);
      m_trainedModelConfigFile[selectedAlgorithm][i] = pathToArch + "/" + filesInSubDir[i];

      // metadata
      auto config = YAML::LoadFile(m_trainedModelConfigFile[selectedAlgorithm][i]);

      auto subjects = new QLabel(config["subjects"].as< std::string >().c_str());
      auto timestamp = new QLabel(config["timestamp"].as< std::string >().c_str());
      auto architecture = new QLabel(config["architecture"].as< std::string >().c_str());
      auto filename = new QLabel(filesInSubDir[i].c_str());

      // put things in order
      trainedModelsLayout->addWidget(m_algoNestedCheckBoxMap[selectedAlgorithm][i], i + 1, 0, Qt::AlignCenter);
      trainedModelsLayout->addWidget(subjects, i + 1, 1, Qt::AlignCenter);
      trainedModelsLayout->addWidget(timestamp, i + 1, 2, Qt::AlignCenter);
      trainedModelsLayout->addWidget(architecture, i + 1, 3, Qt::AlignCenter);
      trainedModelsLayout->addWidget(filename, i + 1, 4, Qt::AlignCenter);
    }

    auto buttonWidth = QLabel().fontMetrics().width("ButtonSize------------------");
    auto buttonHeight = QLabel().fontMetrics().width("----------");

    // initialize the ok and cancel buttons and make them do something meaningful
    auto okButton = new QPushButton("Ok");
    auto cancelButton = new QPushButton("Cancel");

    okButton->setFixedWidth(buttonHeight * 2);
    cancelButton->setFixedWidth(buttonHeight * 2);
    trainedModelsLayout->addWidget(okButton, filesInSubDir.size() + 1, 3, Qt::AlignCenter);
    trainedModelsLayout->addWidget(cancelButton, filesInSubDir.size() + 1, 4, Qt::AlignCenter);
    connect(okButton, SIGNAL(clicked()), this, SLOT(onOkButtonClicked()));
    connect(cancelButton, SIGNAL(clicked()), this, SLOT(onCancelButtonClicked()));
    
    m_dialog->setMaximumSize(buttonWidth * 4, buttonHeight * filesInSubDir.size());
    m_dialog->setLayout(trainedModelsLayout);
    m_dialog->setWindowModality(Qt::ApplicationModal);
    m_dialog->exec();
  }
}

void fSegmentationPanel::onOkButtonClicked()
{
  int counter = 0;
  for (size_t i = 0; i < m_algoNestedCheckBoxMap[m_currentSegmentationAlgorithm].size(); i++)
  {
    if (m_algoNestedCheckBoxMap[m_currentSegmentationAlgorithm][i]->isChecked())
    {
      counter++;
    }
  }
  m_algoCounterMap[m_currentSegmentationAlgorithm]->setText(QString::number(counter));
  // no need to do anything here; all control is happening in getSelectedSegmentationAlgorithms()
  m_dialog->close();
}

void fSegmentationPanel::onCancelButtonClicked()
{
  for (size_t i = 0; i < m_algoNestedCheckBoxMap[m_currentSegmentationAlgorithm].size(); i++)
  {
    m_algoNestedCheckBoxMap[m_currentSegmentationAlgorithm][i]->setChecked(false);
  }
  m_algoCheckBoxMap[m_currentSegmentationAlgorithm]->setChecked(false);
  m_dialog->close();
}

void fSegmentationPanel::onComputeButtonClicked()
{
  emit m_btnComputeClicked();
}

void fSegmentationPanel::helpClicked()
{
  emit helpClicked_SegmentationUsage("gs_drawing.html");
}

std::vector< std::string > fSegmentationPanel::getSelectedSegmentationAlgorithms()
{
  std::vector< std::string > returnVector;
  for (size_t i = 0; i < m_algoCheckBoxMap.size(); i++) // loop through all the segmentation algorithms
  {
    if (m_algoCheckBoxMap[i]->isChecked())
    {
      for (size_t j = 0; j < m_algoNestedCheckBoxMap[i].size(); j++)
      {
        if (m_algoNestedCheckBoxMap[i][j]->isChecked())
        {
          returnVector.push_back(m_trainedModelConfigFile[i][j]);
        }
      }
    }
  }

  return returnVector;
}

std::vector< std::string > fSegmentationPanel::getSelectedFusionAlgorithms()
{
  std::vector< std::string > returnVector;
  for (auto it = m_fusionCheckBoxMap.begin(); it != m_fusionCheckBoxMap.end(); ++it)
  {
    if (it->second->isChecked())
    {
      returnVector.push_back(CAPTK::FusionAlgorithmsString[it->first]);
    }
  }
  return returnVector;
}