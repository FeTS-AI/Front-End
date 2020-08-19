
#ifndef UI_fSegmentationPanel_H
#define UI_fSegmentationPanel_H

#include <QtCore/QVariant>
// #include <QtGui/QAction>
// #include <QtGui/QApplication>
// #include <QtGui/QButtonGroup>
// #include <QtGui/QFrame>
// #include <QtGui/QGridLayout>
// #include <QtGui/QHBoxLayout>
// #include <QtGui/QHeaderView>
// #include <QtGui/QLabel>
// #include <QtGui/QPushButton>
// #include <QtGui/QRadioButton>
// #include <QtGui/QSpacerItem>
// #include <QtGui/QTableWidget>
// #include <QtGui/QVBoxLayout>
// #include <QtGui/QWidget>
// NEW CHANGES
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QFrame>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QRadioButton>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QTableWidget>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>
#include <qcombobox.h>
#include <QToolButton>
#include <qgroupbox.h>
#include <qsize.h>
#include <QCheckBox>
#include <QLineEdit>
#include <map>
//#include "CAPTk.h"
#include "CaPTkGUIUtils.h"
#include "CaPTkEnums.h"

/**
\class QCustomCheckBox

\brief A custom CheckBox class that only sends signal if active
*/
class QCustomCheckBox : public QCheckBox
{
  Q_OBJECT
public:

  QCustomCheckBox(const QString &text, QWidget *parent = 0)
    : QCheckBox(text, parent)
  {
    // Connect original signal to slot 
    connect(this, SIGNAL(clicked()), this, SLOT(detectClicked()));
  }

  ~QCustomCheckBox() {};

signals:
  void clickedOnActive();
  void clickedOnInActive();

private slots:
  void detectClicked()
  {
    if (this->isChecked())
    {
      emit clickedOnActive();
    }
    else
    {
      emit clickedOnInActive();
    }
  }
};

QT_BEGIN_NAMESPACE

class Ui_fSegmentationPanel
{
public:

  QCheckBox *m_deviceControl;
  QPushButton *HelpButton;
  QPushButton * m_btnCompute;
  QPushButton * m_btnInput;
  QLineEdit* m_txtSaveFileName;
  QLineEdit* m_txtInputDirName;
  
  //QPushButton* m_btnBrowseSaveFile;

  QComboBox *m_availableModels;

  std::string m_tempFolderLocation;
  std::map< int, std::map< int, QCheckBox* > > m_algoNestedCheckBoxMap; //! variable that controls which trained models are selected for inference
  std::map< int, QCustomCheckBox* > m_algoCheckBoxMap; //! variable that controls which architectures are enabled
  std::map< int, QLabel* > m_algoCounterMap; //! variable that shows how many architectures are selected
  std::map< int, QCustomCheckBox* > m_fusionCheckBoxMap; //! variable that controls which fusion methods are enabled

  QDialog* m_dialog; //! dialog for individual popups for segmentations

  QGridLayout* m_segAlgoLayout;

  void setupUi(QWidget *parent)
  {
    // get acceptable size of a button
    int buttonWidth = QLabel().fontMetrics().width("ButtonSize------------------");

    // files in the fets data directory
    auto archsInDir = cbica::subdirectoriesInDirectory(m_fetsDataDir.toStdString(), false, true);

    // Browse input fileName
    QGroupBox* inputGroup = new QGroupBox("Input Directory");
    QGroupBox* inputGroup_file = new QGroupBox("Location");
    QHBoxLayout* input_flLayout = new QHBoxLayout();
    //  QLabel* label = new QLabel("FileName:");

    m_btnInput = new QPushButton("...");
    m_btnInput->setToolTip("Browse to select file");

    m_txtInputDirName = new QLineEdit(std::string(loggerFolder).c_str());
    m_txtInputDirName->setAlignment(Qt::AlignCenter | Qt::AlignVCenter);

    input_flLayout->addWidget(m_txtInputDirName);
    input_flLayout->addWidget(m_btnInput);
    inputGroup_file->setLayout(input_flLayout);

    QVBoxLayout* inputLayout = new QVBoxLayout();
    inputLayout->addWidget(inputGroup_file);
    inputLayout ->addStretch();
    inputGroup->setLayout(inputLayout);
    inputGroup->setMaximumWidth(buttonWidth * 3);

    // layout controls
    QGroupBox* algoGroup = new QGroupBox("Algorithms");
    QHBoxLayout* algoLayout = new QHBoxLayout();
    //QVBoxLayout* segAlgoLayout = new QVBoxLayout();
    QVBoxLayout* fusionAlgoLayout = new QVBoxLayout();

    m_segAlgoLayout = new QGridLayout();
    //auto counterLabel = new QLabel("<b>#</b>");
    //m_segAlgoLayout->addWidget(counterLabel, 0, 1, Qt::AlignCenter);

    if (!archsInDir.empty())
    {
      // populate the segmentation algorithms checkboxes
      for (size_t i = 0; i < archsInDir.size(); i++)
      {
        m_algoCheckBoxMap[i] = new QCustomCheckBox(cbica::getFilenameBase(archsInDir[i], false).c_str());
        m_algoCheckBoxMap[i]->setChecked(false);

        // todo: delete this when more algorithms are integrated
        if (m_algoCheckBoxMap[i]->text() == "3DResUNet")
        {
          m_algoCheckBoxMap[i]->setChecked(true);
        }
        m_algoCheckBoxMap[i]->setEnabled(false);
        m_algoCheckBoxMap[i]->setToolTip("Available soon in Phase 2");

        m_algoCounterMap[i] = new QLabel("0");

        m_segAlgoLayout->addWidget(m_algoCheckBoxMap[i], i/* + 1*/, 0, Qt::AlignLeft);
        m_segAlgoLayout->addWidget(m_algoCounterMap[i], i/* + 1*/, 1, Qt::AlignCenter);
      }
    }

    // populate the fusion algorithm checkboxes
    for (size_t i = 0; i < CAPTK::FusionAlgorithmsEnum::FusAlgosMax; i++)
    {
      m_fusionCheckBoxMap[i] = new QCustomCheckBox(CAPTK::FusionAlgorithmsString[i]);
      m_fusionCheckBoxMap[i]->setChecked(false);

      // todo: delete this when more algorithms are integrated
      m_fusionCheckBoxMap[i]->setEnabled(false);
      m_fusionCheckBoxMap[i]->setToolTip("Available soon in Phase 2");

      fusionAlgoLayout->addWidget(m_fusionCheckBoxMap[i]);
    }

    //m_fusionCheckBoxMap[CAPTK::FusionAlgorithmsEnum::GxBoost]->setChecked(false);
    //m_fusionCheckBoxMap[CAPTK::FusionAlgorithmsEnum::STAPLE]->setChecked(false);

    QGroupBox *segmentationGroup = new QGroupBox(("Architectures"));
    QGroupBox *fusionGroup = new QGroupBox(("Fusion"));

    m_availableModels = new QComboBox(parent);
    m_availableModels->insertItem(0, "Skull Stripper");
    m_availableModels->insertItem(1, "Brain Tumor Segmentation");
    fixComboBox(m_availableModels);

    QHBoxLayout* algoLayouthbox = new QHBoxLayout();
    algoLayouthbox->addWidget(m_availableModels);
    algoLayouthbox->addWidget(segmentationGroup);
    algoLayouthbox->addWidget(fusionGroup);
    
    segmentationGroup->setLayout(m_segAlgoLayout);
    fusionGroup->setLayout(fusionAlgoLayout);

    segmentationGroup->setMinimumWidth(buttonWidth * 1.25);
    fusionGroup->setMinimumWidth(buttonWidth);

    algoLayout->addLayout(algoLayouthbox);

    algoGroup->setLayout(algoLayout);
    
    // Browse output fileName
    QGroupBox* saveGroup = new QGroupBox("Output");
    QGroupBox* saveGroup_file = new QGroupBox("Save Location");
    QHBoxLayout* flLayout = new QHBoxLayout();
  //  QLabel* label = new QLabel("FileName:");
    //m_btnBrowseSaveFile = new QPushButton("...");
    //m_btnBrowseSaveFile->setToolTip("Browse to select file");

    m_deviceControl = new QCheckBox("Run computation on GPU");
    m_deviceControl->setToolTip("If disabled, run on CPU");
    m_deviceControl->setChecked(true);

    // this is done solely for the reason for saving everything in the tempDir
    // this folder is deleted after this command to ensure no conflict with m_tempFolderLocation from fMainWindow
    m_txtSaveFileName = new QLineEdit(std::string(loggerFolder).c_str());
    m_txtSaveFileName->setAlignment(Qt::AlignCenter | Qt::AlignVCenter);
    m_txtSaveFileName->setReadOnly(true);
    m_txtSaveFileName->setToolTip("This cannot be edited and is always the same as the input directory");

    m_btnCompute = new QPushButton(parent);
    m_btnCompute->setText(QString("Fusion + Save"));
    m_btnCompute->setToolTip(QString("Fusion and Save Outputs"));
    m_btnCompute->setFixedWidth(buttonWidth + buttonWidth / 5);

   // flLayout->addWidget(label);
    flLayout->addWidget(m_txtSaveFileName);
    //flLayout->addWidget(m_btnBrowseSaveFile);
    //flLayout->addWidget(csv_format);
    //flLayout->addWidget(xml_format);
    saveGroup_file->setLayout(flLayout);

    QVBoxLayout* saveLayout = new QVBoxLayout();
    //saveLayout->addLayout(flLayout);
    saveLayout->addWidget(m_deviceControl);
    saveLayout->addWidget(saveGroup_file);
    saveLayout->addWidget(m_btnCompute);
    saveLayout->addStretch();
    saveGroup->setLayout(saveLayout);
    saveGroup->setMaximumWidth(buttonWidth * 3);


    HelpButton = new QPushButton();
    std::string iconDir = getCaPTkDataDir() + "/icons/";
    HelpButton->setIcon(QIcon((iconDir + "help.png").c_str()));
    HelpButton->setToolTip("Get Help");
    //HelpButton->setIconSize(QSize(30, 30));
    QVBoxLayout *helpLayout = new QVBoxLayout();
    helpLayout->addWidget(HelpButton);
    helpLayout->addStretch();

    QHBoxLayout* subLayout = new QHBoxLayout();
    subLayout->addWidget(inputGroup);
    subLayout->addWidget(algoGroup);
    //subLayout->addWidget(selectionGroup);
    subLayout->addWidget(saveGroup);
    subLayout->addStretch();
    subLayout->addLayout(helpLayout);

    QVBoxLayout* mainLayout = new QVBoxLayout(parent);
    mainLayout->addLayout(subLayout);
    //mainLayout->addStretch();
  }

  //! so that this is done only once
  QString GetFetsDataDir()
  {
    return m_fetsDataDir;
  }

  private:
    QString m_fetsDataDir = cbica::normalizePath(getCaPTkDataDir() + "/fets").c_str();

};

namespace Ui {
  class fSegmentationPanel : public Ui_fSegmentationPanel {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_fSegmentationPanel_H
