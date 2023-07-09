# FeTS: Federated Tumor Segmentation 

The Federated Tumor Segmentation (FeTS) platform, describes an on-going under development open-source toolkit, with a user-friendly graphical user interface (GUI), aiming at:

- bringing pre-trained segmentation models of numerous deep learning algorithms and their fusion, closer to clinical experts and researchers, thereby enabling easy quantification of new radiographic scans and comparative evaluation of new algorithms.
- allowing secure multi-institutional collaborations via federated learning to improve these pre-trained models without sharing patient data, thereby overcoming legal, privacy, and data-ownership challenges.

This project has lead to an easy-to-use potentially-translatable tool enabling easy, fast, objective, repeatable and accurate tumor segmentation, without requiring a computational background by the user, and while facilitating further analysis of tumor radio-phenotypes towards accelerating discovery. 

FeTS is developed and maintained by the [Center for Biomedical Image Computing and Analytics (CBICA)](https://www.cbica.upenn.edu/) at the University of Pennsylvania, in collaboration with [Intel Labs](https://www.intel.com/content/www/us/en/research/overview.html), [Intel AI](https://www.intel.com/ai), and [Intel Internet of Things Group](https://www.intel.com/iot
).

For more details, please visit us at [www.fets.ai](https://www.fets.ai/).

**NOTE**
- If you are looking for the 2022 FeTS study on pre-operative glioblastoma segmentation, please visit the following link: [fets-ai.github.io/FL-Pre](https://fets-ai.github.io/FL-Pre/)
- If you are looking for the 2023 FeTS study on post-operative glioblastoma segmentation, please visit the following link: 
[fets-ai.github.io/FL-PoST](https://fets-ai.github.io/FL-PoST/)

## Citation

If you have found the code and/or its ideas useful, please cite the following articles:

```
@article{fets_study,
  title={Federated learning enables big data for rare cancer boundary detection},
  author={Pati, Sarthak and Baid, Ujjwal and Edwards, Brandon and Sheller, Micah and Wang, Shih-Han and Reina, G Anthony and Foley, Patrick and Gruzdev, Alexey and Karkada, Deepthi and Davatzikos, Christos and others},
  journal={Nature communications},
  volume={13},
  number={1},
  pages={7346},
  year={2022},
  publisher={Nature Publishing Group UK London},
  doi={10.1038/s41467-022-33407-5}
}

@article{fets_tool,
	author={Pati, Sarthak and Baid, Ujjwal and Edwards, Brandon and Sheller, Micah J and Foley, Patrick and Reina, G Anthony and Thakur, Siddhesh P and Sako, Chiharu and Bilello, Michel and Davatzikos, Christos and Martin, Jason and Shah, Prashant and Menze, Bjoern and Bakas, Spyridon},
	title={The federated tumor segmentation (FeTS) tool: an open-source solution to further solid tumor research},
	journal={Physics in Medicine \& Biology},
	url={http://iopscience.iop.org/article/10.1088/1361-6560/ac9449},
	doi={10.1088/1361-6560/ac9449},
	year={2022},
	publisher={IOP Publishing}
	abstract={Objective: De-centralized data analysis becomes an increasingly preferred option in the healthcare domain, as it alleviates the need for sharing primary patient data across collaborating institutions. This highlights the need for consistent harmonized data curation, pre-processing, and identification of regions of interest based on uniform criteria. Approach: Towards this end, this manuscript describes the \textbf{Fe}derated \textbf{T}umor \textbf{S}egmentation (FeTS) tool, in terms of software architecture and functionality. Main Results: The primary aim of the FeTS tool is to facilitate this harmonized processing and the generation of gold standard reference labels for tumor sub-compartments on brain magnetic resonance imaging, and further enable federated training of a tumor sub-compartment delineation model across numerous sites distributed across the globe, without the need to share patient data. Significance: Building upon existing open-source tools such as the Insight Toolkit (ITK) and Qt, the FeTS tool is designed to enable training deep learning models targeting tumor delineation in either centralized or federated settings. The target audience of the FeTS tool is primarily the computational researcher interested in developing federated learning models, and interested in joining a global federation towards this effort. The tool is open sourced at https://github.com/FETS-AI/Front-End.}
}

@software{fets_frontend,
  author       = {Sarthak Pati and
                  Spyridon (Spyros) Bakas},
  title        = {FETS-AI/Front-End: Release for zenodo},
  month        = aug,
  year         = 2022,
  publisher    = {Zenodo},
  version      = {0.0.8},
  doi          = {10.5281/zenodo.7036038},
  url          = {https://doi.org/10.5281/zenodo.7036038}
}
```

## Consensus Models for 2021 Study

Currently, these are available through the FeTS Installer (see Downloads). The models can be downloaded from [this link](https://upenn.box.com/v/fets2021consensusmodels).

## Table of Contents
- [Application Setup](./setup.md)
- [Process the Data](./process_data.md)
- [Run the Application](./runningApplication.md)
- [Extras](./extras.md)
