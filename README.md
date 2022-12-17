# FeTS: Federated Tumor Segmentation 

The Federated Tumor Segmentation (FeTS) platform, describes an on-going under development open-source toolkit, with a user-friendly graphical user interface (GUI), aiming at:

- bringing pre-trained segmentation models of numerous deep learning algorithms and their fusion, closer to clinical experts and researchers, thereby enabling easy quantification of new radiographic scans and comparative evaluation of new algorithms.
- allowing secure multi-institutional collaborations via federated learning to improve these pre-trained models without sharing patient data, thereby overcoming legal, privacy, and data-ownership challenges.

This project has resulted in the development of an easy-to-use translatable tool enabling easy, fast, objective, repeatable and accurate tumor segmentation, without requiring a computational background by the user, and while facilitating further analysis of tumor radio-phenotypes towards accelerating discovery. 

FeTS is developed and maintained by the [Center for Biomedical Image Computing and Analytics (CBICA)](https://www.cbica.upenn.edu/) at the University of Pennsylvania, in collaboration with [Intel Labs](https://www.intel.com/content/www/us/en/research/overview.html), [Intel AI](https://www.intel.com/ai), and Intel Internet of Things Group.

For more details, please visit us at https://www.fets.ai/

## Citation

If you have found the code and/or its ideas useful, please cite the following articles:

```
@article{pati2022federated,
  title={Federated learning enables big data for rare cancer boundary detection},
  author={Pati, Sarthak and Baid, Ujjwal and Edwards, Brandon and Sheller, Micah and Wang, Shih-Han and Reina, G Anthony and Foley, Patrick and Gruzdev, Alexey and Karkada, Deepthi and Davatzikos, Christos and others},
  journal={Nature Communications},
  volume={13},
  number={1},
  pages={1--17},
  year={2022},
  publisher={Nature Publishing Group},
  doi={10.1038/s41467-022-33407-5}
}

@article{pati2022federated,
  title={The federated tumor segmentation (FeTS) tool: an open-source solution to further solid tumor research},
  author={Pati, Sarthak and Baid, Ujjwal and Edwards, Brandon and Sheller, Micah J and Foley, Patrick and Reina, G Anthony and Thakur, Siddhesh and Sako, Chiharu and Bilello, Michel and Davatzikos, Christos and others},
  journal={Physics in Medicine \& Biology},
  volume={67},
  number={20},
  pages={204002},
  year={2022},
  publisher={IOP Publishing},
  doi={10.1088/1361-6560/ac9449}
}

@misc{pati2021gandlf,
      title={GaNDLF: A Generally Nuanced Deep Learning Framework for Scalable End-to-End Clinical Workflows in Medical Imaging}, 
      author={Sarthak Pati and Siddhesh P. Thakur and Megh Bhalerao and Spyridon Thermos and Ujjwal Baid and Karol Gotkowski and Camila Gonzalez and Orhun Guley and Ibrahim Ethem Hamamci and Sezgin Er and Caleb Grenko and Brandon Edwards and Micah Sheller and Jose Agraz and Bhakti Baheti and Vishnu Bashyam and Parth Sharma and Babak Haghighi and Aimilia Gastounioti and Mark Bergman and Anirban Mukhopadhyay and Sotirios A. Tsaftaris and Bjoern Menze and Despina Kontos and Christos Davatzikos and Spyridon Bakas},
      year={2021},
      eprint={2103.01006},
      archivePrefix={arXiv},
      primaryClass={cs.LG}
}
```

## Supporting Grant

This work is in part supported by the National Institutes of Health / National Cancer Institute / Informatics Technology for Cancer Research (NIH/NCI/ITCR), under grant award number U01-CA242871.

## Downloads

By downloading FeTS, you agree to our [License](./LICENSE). 

Latest installer: https://www.nitrc.org/frs/downloadlink.php/13219

## Documentation

https://fets-ai.github.io/Front-End/

## Disclaimer

- The software has been designed for research purposes only and has neither been reviewed nor approved for clinical use by the Food and Drug Administration (FDA) or by any other federal/state agency.
- Certain part of the code for this user interface (excluding dependent libraries) is governed by the license provided in https://www.med.upenn.edu/sbia/software-agreement.html unless otherwise specified.

For more details, please visit us at https://www.fets.ai/

For issues, please visit https://github.com/FETS-AI/Front-End/issues 

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
