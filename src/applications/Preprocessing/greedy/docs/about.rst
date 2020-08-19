***************
What is Greedy?
***************

Greedy is a tool for fast medical image registration. It was developed by Paul Yushkevich at the `Penn Image Computing and Science Lab` at the University of Pennsylvania. The motivation for developing greedy was to have a really fast CPU-based deformable image registration tool that could be used in applications where many images have to be registered in parallel - like multi-atlas image segmentation. 

Greedy shares multiple concepts implementation strategies with `Symmetric Normalization (SyN) in ANTS`. But greedy is non-symmetric, which makes it faster (in applications like multi-atlas segmentation, symmetric property is not required). Greedy also uses highly optimized code for image metric computation that adds extra speed. 
This work is funded by the NIH/NIBIB under grants R01 EB-017255 and R01 EB-014146

.. _Penn Image Computing and Science Lab: http://picsl.upenn.edu/
.. _Symmetric Normalization (SyN) in ANTS: http://stnava.github.io/ANTs/
