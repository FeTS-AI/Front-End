***********
Quick Start
***********

In this example, we will register brain MRI scans from two different individuals. One MRI will be designated as "fixed" (``fixed.nii.gz``) and the other as "moving" (``moving.nii.gz``). We also assume that we have a multi-label segmentation of the moving image, called ``moving_seg.nii.gz``. 

To run the code in this section, you must make sure that the executable ``greedy`` is in your system path. On Linux and MacOS, run::

    > export PATH=/path/to/directory/containing/greedy:$PATH

On Windows, run::

    > PATH=c:\path\to\directory\containing\greedy;%PATH%

Affine Registration
~~~~~~~~~~~~~~~~~~~

The two brains are in different physical locations. We first perform affine registration to correct for differences in brain position, rotation, and size::

    > greedy -d 3 -a \
        -m NCC 2x2x2 \
        -i fixed.nii.gz moving.nii.gz \
        -o affine.mat \
        -ia-image-centers -n 100x50x10

This call to greedy uses some of the most common options:

* ``-d 3`` specifies the dimensionality of the problem (3D registration)
* ``-a`` specifies that we are performing affine registration
* ``-i`` specifies the fixed/moving image pair
* ``-m NCC 2x2x2`` specifies the image dissimilarity metric. Greedy will use the normalized cross-correlation metric with a 2x2x2 patch radius (patch size 5x5x5)
* ``-o`` specifies the file where the output affine transformation (a matrix) will be written
* ``-ia-image-centers`` specifies that the affine transform is initialized by matching the centers of the images. This is useful when your images do not occupy the same physical space.
* ``-n 100x50x10`` instructs greedy to run for 100 iterations at the lowest resolution level, 50 at intermediate resolution and 10 at full resolution.
 

The result of the registration is the file ``affine.mat``. We can view the contents of this file::

    > cat affine.mat

So far, we have not actually applied the registration to the moving image; this requires a separate call to greedy, as discussed below.


Deformable Registration
~~~~~~~~~~~~~~~~~~~~~~~

We can now perform deformable registration between the fixed and moving images::

    > greedy -d 3 \
        -m NCC 2x2x2 \
        -i fixed.nii.gz moving.nii.gz \
        -it affine.mat \
        -o warp.nii.gz \
        -oinv inverse_warp.nii.gz \
        -n 100x50x10

This is pretty similar to the affine command. Notice the differences:

* We are no longer supplying the ``-a`` flag. This activates the default deformable registration mode.

* ``-it`` is used to specify the initial transformation of the moving image. We pass the affine registration output as the initial transformation.

* ``-o`` is used to specify a filename of the image where the resulting dense deformable transformation (warp) will be stored.

* Optional command ``-oinv`` is used to tell greedy to approximate the inverse warp. Inverse warps are rarely needed in practical applications.


After you run the command above, the ``warp.nii.gz`` file will be generated. You can view it in ITK-SNAP or other viewer. It is a 3D image volume where each voxel stores the *x,y,z* components of a displacement vector. *For every voxel in the fixed image space, the warp image specifies the displacement that maps that voxel's center to the corresponding location in the moving image.*

Applying warps to images and segmentations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

So far we have just computed the affine and deformable transformation between the fixed image and the moving image. We now need to apply these transformations to the moving image. We also have a segmentation of the moving image, moving\_seg.nii.gz, which we would also like to warp into the fixed image space. We do this using the greedy's **reslice mode** (``-r`` commands)::

    > greedy -d 3 \
      -rf fixed.nii.gz \
      -rm moving.nii.gz resliced.nii.gz \
      -ri LABEL 0.2vox \
      -rm moving_seg.nii.gz resliced_seg.nii.gz \
      -r warp.nii.gz affine.mat

Let's deconstruct this command:

-  ``-rf`` specifies the referenced/fixed image space. Images will be re-sliced into this space.

-  ``-rm`` specifies an input/output pair. The input is in the moving image space, the output will contain the image resliced (warped) into the fixed image space. Multiple ``-rm`` commands can be provided as above.

-  ``-ri`` specifies the interpolation mode for subsequent ``-rm`` commands. Default interpolation is ``LINEAR`` interpolation. When warping multi-label segmentations, greedy provides a special ``LABEL`` interpolation mode. This mode applies a little bit of smoothing to each label in your segmentation (including the background), warps this smoothed segmentation, and then performs voting among warped smoothed binary segmentations to assign each voxel in reference space a label. This works better than nearest neighbor interpolation (less aliasing). The ``0.2vox`` part of the command specifies the amount of smoothing.

-  ``-r`` is used to specify a sequence of transformations, from last to first.

Applying the Inverse Warp
~~~~~~~~~~~~~~~~~~~~~~~~~

Optionally, you can warp the fixed image back into the space of the moving. This requires using the (``-oinv``) option when calling greedy for deformable registration. You can also invert a warp computed previously using the ``-iw`` command (see documentation).

The following command warps the fixed image to the moving image space::

    > greedy -d 3 \
        -rf moving.nii.gz \
        -rm fixed.nii.gz reslice_fix_into_mov.nii.gz \
        -r affine.mat,-1 inverse_warp.nii.gz

The call to greedy is very similar to above, except that:

-  the roles of fixed and moving is switched (obviously)

-  the affine and deformable transformations are provided in reverse order

-  the affine transformation is inverted (the -1 after ``affine.mat``)

That's pretty much it for learning to use greedy. See detailed documentation for other options and more details. Enjoy!
