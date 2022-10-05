# FeTS: Federated Tumor Segmentation 

The Federated Tumor Segmentation (FeTS) platform, describes an on-going under development open-source toolkit, with a user-friendly graphical user interface (GUI), aiming at:

- bringing pre-trained segmentation models of numerous deep learning algorithms and their fusion, closer to clinical experts and researchers, thereby enabling easy quantification of new radiographic scans and comparative evaluation of new algorithms.
- allowing secure multi-institutional collaborations via federated learning to improve these pre-trained models without sharing patient data, thereby overcoming legal, privacy, and data-ownership challenges.

Successful completion of this project will lead to an easy-to-use potentially-translatable tool enabling easy, fast, objective, repeatable and accurate tumor segmentation, without requiring a computational background by the user, and while facilitating further analysis of tumor radio-phenotypes towards accelerating discovery. 

FeTS is developed and maintained by the [Center for Biomedical Image Computing and Analytics (CBICA)](https://www.cbica.upenn.edu/) at the University of Pennsylvania, in collaboration with [Intel Labs](https://www.intel.com/content/www/us/en/research/overview.html), [Intel AI](https://www.intel.com/ai), and Intel Internet of Things Group.

For more details, please visit us at https://www.fets.ai/

## Citation

If you have found the code and/or its ideas useful, please cite the following articles:

```
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

## Supporting Grant

This work is in part supported by the National Institutes of Health / National Cancer Institute / Informatics Technology for Cancer Research (NIH/NCI/ITCR), under grant award number U01-CA242871.

## Documentation

https://fets-ai.github.io/Front-End/

## Disclaimer

- The software has been designed for research purposes only and has neither been reviewed nor approved for clinical use by the Food and Drug Administration (FDA) or by any other federal/state agency.
- Certain part of the code for this user interface (excluding dependent libraries) is governed by the license provided in https://www.med.upenn.edu/sbia/software-agreement.html unless otherwise specified.

For more details, please visit us at https://www.fets.ai/

For issues, please visit https://github.com/FETS-AI/Front-End/issues 

## Downloads

By downloading FeTS, you agree to our [License](./LICENSE). 

## Contact

For more information, please contact <a href="mailto:software@cbica.upenn.edu">CBICA Software</a>.

## GitHub Distribution

We currently provide only our tagged versions of the code via GitHub. Check the "tags" using your favorite Git client after cloning our repository. The analogous commands are as follows:

```bash
git clone https://github.com/FETS-AI/Front-End.git
latesttag=$(git describe --tags)
echo checking out ${latesttag}
git checkout ${latesttag}
```
