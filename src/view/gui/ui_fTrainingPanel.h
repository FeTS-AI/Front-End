
#ifndef UI_fTrainingPanel_H
#define UI_fTrainingPanel_H

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
\class QCustomRadioButton

\brief A custom PushButton class that only sends signal if active
*/
class QCustomRadioButton : public QRadioButton
{
Q_OBJECT
public:
  QCustomRadioButton(const QString &text, QWidget *parent = 0)
    : QRadioButton(text, parent)
  {
    // Connect original signal to slot 
    connect(this, SIGNAL(toggled(bool)), this, SLOT(detectClicked()));
  }

  ~QCustomRadioButton() {};

signals:
  void toggledOnActive();

private slots:
  void detectClicked()
  {
    if (this->isChecked())
    {
      emit toggledOnActive();
    }
  }
};

QT_BEGIN_NAMESPACE

class Ui_fTrainingPanel
{
public:

  QCheckBox *m_deviceControl;
  QPushButton *HelpButton;
  QPushButton * m_btnCompute;
  QPushButton * m_btnInput;
  QLineEdit* m_txtSaveFileName;
  QLineEdit* m_txtInputDirName;
  QLineEdit* m_collaboratorName;
  
  //QPushButton* m_btnBrowseSaveFile;

  QDialog* m_dialog; //! dialog for individual popups for segmentations

  QComboBox *m_availableModels;

  std::string m_tempFolderLocation;
  std::map< int, std::map< int, QRadioButton* > > m_algoNestedRadioButtonMap; //! variable that controls which trained models are selected for inference
  std::map< int, QCustomRadioButton* > m_algoRadioButtonMap; //! variable that controls which architectures are enabled

  void setupUi(QWidget *parent)
  {
    int buttonWidth = QLabel().fontMetrics().width("ButtonSize------------------");

    // Browse input fileName
    QGroupBox* inputGroup = new QGroupBox("Inputs");
    QGroupBox* inputGroup_file = new QGroupBox("Data Directory");
    QHBoxLayout* input_flLayout = new QHBoxLayout();
    //  QLabel* label = new QLabel("FileName:");
    m_btnInput = new QPushButton("...");
    m_btnInput->setToolTip("Browse to select file");

    m_txtInputDirName = new QLineEdit(std::string(loggerFolder).c_str());
    m_txtInputDirName->setAlignment(Qt::AlignCenter | Qt::AlignVCenter);

    input_flLayout->addWidget(m_txtInputDirName);
    input_flLayout->addWidget(m_btnInput);
    inputGroup_file->setLayout(input_flLayout);

    // common collaborator name
    QGroupBox* inputGroup_colName = new QGroupBox("Common Collaborator Name");
    QHBoxLayout* input_flLayout_name = new QHBoxLayout();
    //  QLabel* label = new QLabel("FileName:");

    m_collaboratorName = new QLineEdit("Insert the name as it appears in certificate");
    m_collaboratorName->setAlignment(Qt::AlignCenter | Qt::AlignVCenter);

    input_flLayout_name->addWidget(m_collaboratorName);
    inputGroup_colName->setLayout(input_flLayout_name);

    QVBoxLayout* inputLayout = new QVBoxLayout();
    inputLayout->addWidget(inputGroup_file);
    inputLayout->addWidget(inputGroup_colName);
    inputLayout->addStretch();
    inputGroup->setLayout(inputLayout);
    inputGroup->setMaximumWidth(buttonWidth * 3);

    QGroupBox* algoGroup = new QGroupBox("Algorithms");
    QVBoxLayout* algoLayout = new QVBoxLayout();

    m_availableModels = new QComboBox(parent);
    m_availableModels->insertItem(0, "Brain Tumor Segmentation");
    fixComboBox(m_availableModels);
    algoLayout->addWidget(m_availableModels);

    for (size_t i = 0; i < CAPTK::TrainingSegmentationAlgorithmsEnum::train_SegAlgosMax; i++)
    {
      m_algoRadioButtonMap[i] = new QCustomRadioButton(CAPTK::TrainingSegmentationAlgorithmsString[i]);
      // todo: delete this when more algorithms are integrated
      if (m_algoRadioButtonMap[i]->text() == "3DResUNet")
      {
        m_algoRadioButtonMap[i]->setChecked(true);
      }
      m_algoRadioButtonMap[i]->setEnabled(false);
      m_algoRadioButtonMap[i]->setToolTip("Available soon in Phase 2");

      algoLayout->addWidget(m_algoRadioButtonMap[i]);
    }

    algoGroup->setLayout(algoLayout);
    algoGroup->setMinimumWidth(buttonWidth);

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
    m_btnCompute->setText(QString("Train + Save"));
    m_btnCompute->setToolTip(QString("Train and Save Model"));
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

  private:
    QString m_fetsDataDir = cbica::normalizePath(getCaPTkDataDir() + "/fets").c_str();

};

namespace Ui {
  class fTrainingPanel : public Ui_fTrainingPanel {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_fSegmentationPanel_H
