*********
Reference
*********

Greedy Usage
------------

.. code-block:: bash

      Usage:
        greedy [options]
      Required options:
        -d DIM                 : Number of image dimensions
        -i fix.nii mov.nii     : Image pair (may be repeated)
        -o output.nii          : Output file
      Mode specification:
        -a                     : Perform affine registration and save to output (-o)
        -brute radius          : Perform a brute force search around each voxel
        -moments <1|2>         : Perform moments of inertia rigid alignment of given order.
                                     order 1 matches center of mass only
                                     order 2 matches second-order moments of inertia tensors
        -r [tran_spec]         : Reslice images instead of doing registration
                                     tran_spec is a series of warps, affine matrices
        -iw inwarp outwarp     : Invert previously computed warp
        -root inwarp outwarp N : Convert 2^N-th root of a warp
        -jac inwarp outjac     : Compute the Jacobian determinant of the warp
      Options in deformable / affine mode:
        -w weight              : weight of the next -i pair
        -m metric              : metric for the entire registration
                                     SSD:          sum of square differences (default)
                                     MI:           mutual information
                                     NMI:          normalized mutual information
                                     NCC <radius>: normalized cross-correlation
                                     MAHAL:        Mahalanobis distance to target warp
        -e epsilon             : step size (default = 1.0),
                                     may also be specified per level (e.g. 0.3x0.1)
        -n NxNxN               : number of iterations per level of multi-res (100x100)
        -threads N             : set the number of allowed concurrent threads
        -gm mask.nii           : mask for gradient computation
        -gm-trim <radius>      : generate mask for gradient computation by trimming the extent
                                 of the fixed image by given radius. This is useful during affine
                                 registration with the NCC metric when the background of your images
                                 is non-zero. The radius should match that of the NCC metric.  -mm mask.nii           : mask for the moving image
        -it filenames          : sequence of transforms to apply to the moving image first
      Specific to deformable mode:
        -tscale MODE           : time step behavior mode: CONST, SCALE [def], SCALEDOWN
        -s sigma1 sigma2       : smoothing for the greedy update step. Must specify units,
                                 either `vox` or `mm`. Default: 1.732vox, 0.7071vox
        -oinv image.nii        : compute and write the inverse of the warp field into image.nii
        -oroot image.nii       : compute and write the (2^N-th) root of the warp field into image.nii, where
                                 N is the value of the -exp option. In stational velocity mode, it is advised
                                 to output the root warp, since it is used internally to represent the deformation
        -wp VALUE              : Saved warp precision (in voxels; def=0.1; 0 for no compression).
        -noise VALUE           : Standard deviation of white noise added to moving/fixed images when
                                 using NCC metric. Relative to intensity range. Def=0.001
        -exp N                 : The exponent used for warp inversion, root computation, and in stationary
                                 velocity field (Diff Demons) mode. N is a positive integer (default = 6)
        -sv                    : Performs registration using the stationary velocity model, similar to diffeomoprhic
                                 Demons (Vercauteren 2008 MICCAI). Internally, the deformation field is
                                 represented as 2^N self-compositions of a small deformation and
                                 greedy updates are applied to this deformation. N is specified with the -exp
                                 option (6 is a good number). This mode results in better behaved
                                 deformation fields and Jacobians than the pure greedy approach.
        -svlb                  : Same as -sv but uses the more accurate but also more expensive
                                 update of v, v <- v + u + [v,u]. Experimental feature
        -id image.nii          : Specifies the initial warp to start iteration from. In stationary mode, this
                                 is the initial stationary velocity field (output by -oroot option)
      Initial transform specification:
        -ia filename           : initial affine matrix for optimization (not the same as -it)
        -ia-identity           : initialize affine matrix based on NIFTI headers
        -ia-image-centers      : initialize affine matrix based on matching image centers
        -ia-image-side CODE    : initialize affine matrix based on matching center of one image side
        -ia-moments <1|2>      : initialize affine matrix based on matching moments of inertia
      Specific to affine mode (-a):
        -dof N                 : Degrees of freedom for affine reg. 6=rigid, 12=affine
        -jitter sigma          : Jitter (in voxel units) applied to sample points (def: 0.5)
        -search N sa sx [flip] : Random search over rigid transforms (N iter) before starting optimization
                                 sa, sx: sigmas for rot-n angle (degrees) and offset between image centers
                                 'flip' optional keyword will enable search over flips too
      Specific to moments of inertia mode (-moments 2):
        -det <-1|1>            : Force the determinant of transform to be either 1 (no flip) or -1 (flip)
        -cov-id                : Assume identity covariance (match centers and do flips only, no rotation)
      Specific to reslice mode (-r):
        -rf fixed.nii          : fixed image for reslicing
        -rm mov.nii out.nii    : moving/output image pair (may be repeated)
        -rs mov.vtk out.vtk    : moving/output surface pair (vertices are warped from fixed space to moving)
        -ri interp_mode        : interpolation for the next pair (NN, LINEAR*, LABEL sigma)
        -rb value              : background (i.e. outside) intensity for the next pair (default 0)
        -rc outwarp            : write composed transforms to outwarp
        -rj outjacobian        : write Jacobian determinant image to outjacobian
      For developers:
        -debug-deriv           : enable periodic checks of derivatives (debug)
        -debug-deriv-eps       : epsilon for derivative debugging
        -debug-aff-obj         : plot affine objective in neighborhood of -ia matrix
        -dump-moving           : dump moving image at each iter
        -dump-freq N           : dump frequency
        -powell                : use Powell's method instead of LGBFS
        -float                 : use single precision floating point (off by default)
        -version               : print version info

General Options
---------------

Image dimensionality (``-d``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Format: ``-d <2|3>``

This option is required to run greedy. It specifies whether registration
or other operations should be performed in 2D or 3D.

Number of parallel threads (``-threads``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Format: ``-threads <number>``

By default, greedy will run in multithreaded mode, using all of your available CPU cores. You can restrict the number of cores used to any given number. On many clusters, the ``NSLOTS`` environment variable is defined and can be used to set the number of threads correctly::

    > if [[ $NSLOTS -gt 1 ]]; then \
        greedy -d 3 -threads $NSLOTS ... ; \
      fi

Floating point precision (``-float``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Format: ``-float``

By default, greedy uses double precision floating point to represent images and transformations in memory. This option uses single-precision instead. This is faster and uses less memory, but at some small loss of precision (especially during NCC metric computation). *We recommend not using this option, as double precision floating point has been tested far more extensively*.

Command-line help (``-h``)
~~~~~~~~~~~~~~~~~~~~~~~~~~

Format: ``-h``

Use this command to list all the commands and options for greedy. Some commands are esoteric or developer-oriented and are not discussed here.

Common Commands in Deformable Registration Mode
-----------------------------------------------

**Input Image Pair and Weight Specification (-i, -w)**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**-i <fixed\_image> <moving\_image>**

**-w <weight>**

This command specifies the fixed/moving image pair. Multiple such
commands can be provided, in which case there will be multiple fixed and
multiple moving images. However, all the fixed images must be in the
same physical space, as must be all the moving images. You can use the
**-w** command to assign different weights to different fixed/moving
pairs. Note that the **-w** command applies to all subsequent **-i**
commands.

> greedy -d 3 \\

-w 0.25 -i fixed\_t1.nii moving\_t1.nii \\

-w 0.75 -i fixed\_t2.nii moving\_t2.nii \\

...

The fixed and moving images may also be multi-component images (e.g,
images of vectors or tensors).

**Output Warp Specification: (-o)**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**-o <warp\_image>**

Specifies where the warp image will be stored. The warp image will be in
the same space as the fixed image and will have three components per
pixel. The warp image is specified as follows. Suppose that A is a voxel
coordinate in the fixed image and B is a voxel coordinate in the moving
image, and that registration matched A to B. Then

ras(A) + warp[A] = ras(B)

where ras(A) is the physical coordinate of voxel coordinate A in the RAS
coordinate space (space used by NIFTI).

**Metric Specification (-m)**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**-m <SSD \| NMI \| NCC <radius> >**

Specifies the image similarity metric used for the registration. Greedy
does not allow mixing multiple metrics in the same registration
(weighting multiple metrics in non-trivial anyway). So the position of
the command on the command line does not matter.

Three metrics are supported:

-  Sum of squared differences (SSD) - fastest but only suitable for
       same-modality registration where intensity ranges of the fixed
       and moving images are the same. For example, two CT scans. This
       metric just tries to match the intensity of the fixed and moving
       images at every voxel.

-  Normalized cross-correlation (NCC) - relatively fast too, but more
       robust to noise and intensity differences. Tries to maximize the
       correlation coefficient between the neighborhood of each voxel in
       the fixed image and the corresponding neighborhood in the moving
       image. The size of the neighborhood is specified by **<radius>**.
       For example **NCC 2x2x2** specifies a 5x5x5 neighborhood. Note
       that there is almost to performance cost for using larger radii
       due to efficient implementation.

-  Normalized mutual information (NMI) - should be used when intensity
       spaces of the moving and fixed images are very different, e.g.,
       registering T1-MRI to T2-MRI. Does not work very well for
       deformable registration, better for affine/rigid.

> greedy -d 3 \\

-m NCC 4x4x4 -i fixed\_t1.nii moving\_t1.nii \\

...

**Initial Transformations (-it)**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**-it <transform> [transform] ...**

Provides a chain of transformations (affine matrices, warps) to apply to
the moving image before registration. This is equivalent to first
reslicing the moving image into the fixed image space using the same
chain of transformations (**-r** command). The most common scenario is
to provide the output of affine/rigid registration to the **-it**
command.

> greedy -d 3 \\

-it affine.mat -i fixed\_t1.nii moving\_t1.nii \\

...

**Fixed Image Mask (-gm)**
~~~~~~~~~~~~~~~~~~~~~~~~~~

**-gm <mask\_image> **

Specifies a mask that restricts registration to a region of the fixed
image. This can make registration faster and more robust and is highly
recommended, particularly when there is a lot of intensity variation
along the boundaries of the fixed image. The mask image is typically a
binary image, but a soft mask can also be provided.

**Multi-resolution schedule (-n)**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**-n <iteration\_spec> **

Specify how many iterations of registration to do at each iteration
level. For example **-n 100x40x20** does three levels of
super-resolution (4x, 2x and 1x) and does 100 iterations at 4x (coarsest
level), 40 iterations at 2x (intermediate) and 20 iterations at 1x (full
resolution).

**Inverse warp output (-oinv, -invexp)**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**-oinv <warp\_image>**

**-invexp <exponent> **

Unlike symmetric normalization (SyN), greedy does not compute the
inverse of the deformation field at each iteration of image
registration. However, you can still generate an inverse warp post-hoc.
This uses the fixed point method of warp inversion [reference!]. This
adds some extra time at the end of the registration.

To improve the performance of the inverse algorithm, the forward warp is
first taken to a power -2, -4, -8, etc. In other words, we find a warp
psi, such that psi(psi( ... psi(psi(x)))) = warp(x). The **exponent**
parameter to **-invexp** is used to specify the power, with power =
2^-exponent. Default value is 2. If you get bad (self-intersecting)
inverse warps, try a larger value.

> greedy -d 3 \\

... -invexp 4 -oinv inverse\_warp.nii.gz

**Deformable Registration Parameters**
--------------------------------------

**Smoothing Kernels (-s)**
~~~~~~~~~~~~~~~~~~~~~~~~~~

**-s <gradient\_sigma> <warp\_sigma>**

Probably the most crucial parameter for deformable registration. This
specifies the amount of regularization applied to the deformation field
during registration. Just like in SyN (and in Demons registration before
that), there are two types of regularization applied:

-  Metric gradient regularization: this is applied to smooth the
       gradient of the image match metric at each iteration. The
       smoothed gradient is used to update the current estimate of the
       warp via composition. Larger values of smoothing
       (**gradient\_sigma**) result is smoother deformation fields.

   -  The default value of **gradient\_sigma** is 1.732vox (square root
          of 3). This default matches the default in SyN.

-  Warp regularization: the entire warp field is smoothed after each
       iteration. This dampens the overall deformation. Larger values of
       **warp\_sigma** give smaller deformations.

   -  The default value of **gradient\_sigma** is 0.707vox (square root
          of 0.5). This default matches the default in SyN.

Both sigmas can be provided in voxel units (suffix **vox**) or physical
units (suffix **mm**).

> greedy -d 3 -s 2mm 0.7mm ...

> greedy -d 3 -s 1.5x1.8x2.0vox 0.2vox ...

**Step Size (-e) and time step scaling (-tscale)**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**-e <step\_spec>**

**-tscale <SCALE \| SCALEDOWN \| CONST>**

Command **-e** specifies the "time step" size used to update the warp at
each iteration. Larger values can speed up registration but can also
cause deformation to become non-diffeomorphic. The default value is 1.0,
and typical values are in the 0.25 to 0.5 range.

> greedy -d 3 -e 0.5 ...

> greedy -d 3 -n 100x40x20 -e 1.0x0.5x0.2 ...

The second form of the command specifies different step size for each
multi-resolution level. This has not proven useful in my experience.

By default, the time step is applied after scaling the smoothed metric
gradient so that the norm of the largest gradient across the whole image
is 1 voxel. This behavior can be modified with the **-tscale** command,
but this is not recommended. Other options are **SCALEDOWN** (where the
gradient is only scaled down to have maximum norm 1 but never up) and
**CONST** (the gradient is never scaled, so you have to set your time
step extremely carefully).

**Warp field precision (-wp)**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**-e <real\_value>**

Warp fields have great potential to take over disk space. By default,
greedy stores warp fields only to the precision of 1/10 of voxel size.
In most applications, there is no real difference to warping an image by
2.2 voxels or 2.24 voxels. By lowering precision, you can achieve much
better compression when storing warp files in **.nii.gz** and other
compressed formats. You can change the precision from the default 0.1
(1/10 voxel) to full precision (0) or any other value between 0 and 1.

> greedy -d 3 -wp 0.01 ...

**Affine and Rigid Registration**
---------------------------------

**Affine mode (-a, -dof)**
~~~~~~~~~~~~~~~~~~~~~~~~~~

**-a**

**-o <affine\_matrix>**

**-dof <6\|7\|12>**

Calling greedy with **-a** command switches the tool to affine/rigid
mode. Affine/rigid mode can not be combined with deformable mode in the
same command.

By default, full affine registration is performed (12 degrees of freedom
in 3D). To use rigid registration, pass in **-dof 6**. To use rigid +
uniform scaling, use **-dof 7**.

In affine mode, many of the same options as in deformable mode are used,
with some minor differences.

-  **-o** command will write out a matrix encoding the affine transform.
       This is a N+1 x N+1 matrix that maps voxels in fixed image space
       to voxels in moving image space. Specifically, if voxel
       coordinate A in the fixed image corresponds to voxel coordinate B
       in the moving image, then

[ras(B); 1] = Matrix \* [ras(B); 1]

-  If you wish to convert the matrix file to a different format or
       perform various operations on matrix files, use the
       **c3d\_affine\_tool** in
       `*Convert3D* <http://itksnap.org/c3d>`__.

-  **-i**, **-w**, **-m,** **-n, -gm** behave the same way as in defor
       mable mode.

-  **-ia** or **-ia-identity** should be used to initialize affine
       registration (instead of **-it** in deformable mode)

-  **-s** and **-e** have no effect.

-  **-oinv** is not supported. If you want to invert the affine
       transformation, use the **c3d\_affine\_tool** in
       `*Convert3D* <http://itksnap.org/c3d>`__.

Typical example of rigid registration:

> greedy -d 3 **-a** \\

-i fixed.nii.gz moving.nii.gz \\

-gm fixed\_mask.nii.gz \\

-ia-identity \\

-dof 6 -o rigid.mat \\

-n 100x50x0 -m NCC 4x4x4

**Initial transform specification for affine/rigid mode**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**-ia <affine\_matrix>**

**-ia-identity**

**-ia-image-centers**

You can initialize rigid/affine registration with a given matrix or with
the identity matrix. Using the identity matrix will initialize the image
alignment based on image headers (i.e., assume that ras(A) = ras(B)).
Command **-ia-image-centers** matches image centers (by translation).

If **neither** of these three options is given, images are initialized
based on voxel coordinates, rather than on image headers. This can
result in registration failures for many images.

**Affine initialization via random search (-search)**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**-search <N\_iter> <sigma\_angle> <sigma\_offset>**

This command will randomly sample **N\_iter** starting positions for
affine registration and start optimization from the best position found.
Random sampling generates rigid transformations of the moving image. The
**sigma** parameters specify the range of the angles of rotation (in
degrees) and range of the offset (in voxels).

> greedy -d 3 -a \\

-i fixed.nii.gz moving.nii.gz \\

-gm fixed\_mask.nii.gz \\

-ia-identity \\

-dof 6 -o rigid.mat \\

-n 100x50x0 -m NCC 4x4x4 \\

**-search 1000 10 20**

**Random jitter (-jitter)**
~~~~~~~~~~~~~~~~~~~~~~~~~~~

**-jitter <real\_value>**

Affine registration tends to converge better when the sample locations
where metric is calculated are randomly displaced from voxel centers
(this avoids spurious local minima). By default a random jitter with
range [-0.5 0.5] is applied to the voxel coordinates where images are
sampled. For faster initialization, set jitter to 0.0.

**Image Reslicing Mode**
------------------------

The image reslicing mode is used to apply warps and affine matrices to
images. It can also be used to compose multiple transforms into a single
transform, and to apply warps to meshes. Reslicing mode is activated
when the **-r** command is used. Reslicing mode cannot be combined with
registration in the same command line.

-  See examples under *Quick Start*

**Reference (fixed) space specification (-rf)**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**-rf <reference\_image>**

Specify the reference image for the reslicing. All images will be
resliced into the space of the reference image.

**Input/output pair specification (-rm)**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**-rm <input\_image> <output\_image>**

Specify an image to be resliced and the corresponding output image. You
can have any number of **-rm** commands in the same command line. The
input images provided to **-rm** commands do not have to be in the same
physical space. They can be scalar or multi-component images.

**Interpolation mode specification (-ri)**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**-ri <NN \| LINEAR \| LABEL <sigma\_spec> >**

Specify the interpolation mode to use for reslicing. This command only
affects the subsequent **-rm** commands on the command line (so should
precede the **-rm** command you want it to affect). The following modes
are available;

-  Nearest neighbor (NN): rarely recommended, results in the most
       aliasing

-  Bilinear/trilinear interpolation (LINEAR): default interpolation
       mode, fast and less aliasing

-  Label interpolation (LABEL): this special mode is used for
       warping/reslicing multi-label segmentations. This mode applies a
       little bit of smoothing to each label in your segmentation
       (including the background), warps this smoothed segmentation, and
       then performs voting among warped smoothed binary segmentations
       to assign each voxel in reference space a label. This works
       better than nearest neighbor interpolation (less aliasing).

   -  The **<sigma\_spec>** parameter to the **-ri LABEL** command
          specifies the standard deviation of the Gaussian kernel used
          to smooth the labels. It can be provided in voxel units (e.g.,
          0.2vox) or millimeter units (e.g., 0.2mm). Value of 0.2vox
          works well in most situations.

> greedy -d 3 \\

-rf reference.nii \\

-ri LINEAR \\

-rm t1mri.nii.gz warped\_t1mri.nii.gz \\

-ri LABEL 0.2vox \\

-rm segmentation.nii.gz warped\_seg.nii.gz \\

-r ...

**Transform chain specification (-r)**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**-r <transform\_spec> [transform\_spec] ...**

Specify the chain of transforms to be applied to the moving image. The
last transform is applied first. In most cases, you would do affine
registration followed by deformable registration. To reslice your
original moving image into the space of the fixed image you would use
the command

> greedy -d 3 \\

-rf fixed.nii.gz \\

-rm moving.nii.gz resliced.nii.gz \\

-r warp.nii.gz affine.mat

So the moving image will first be transformed by the affine transform,
and then by the warp. Or in other words, if A is a voxel coordinate in
fixed image space, then the corresponding voxel coordinate B in the
moving image is found according to

ras[B] = warp(affine(ras[A])

-  It is a common error to provide transforms in the wrong order.

-  You can provide as many transforms as you wish - it is possible to
       chain a dozen transforms.

-  To specify that the affine transform should be inverted, use
       **affine.mat,-1** syntax.

   -  For example, to reslice the fixed image into the space of the
          moving image in the above example, use

> greedy -d 3 \\

-rf moving.nii.gz \\

-rm fixed.nii.gz resliced\_backwards.nii.gz \\

-r affine.mat,-1 inverse\_warp.nii.gz

-  Note that the order of transforms has switched. This is because

ras[A] = inverse\_affine(inverse\_warp((ras[B])

**Composing transformations (-rc)**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**-rc <warp\_image>**

In addition (or instead of) reslicing images, you can use the reslice
mode to compose multiple transforms or to convert an affine transform
into the corresponing (linear) warp image. For example:

> greedy -d 3 \\

-rf fixed.nii.gz \\

-r warp1.nii.gz warp2.nii.gz affine.mat \\

-rc composite\_warp.nii.gz

The **-rc** command can be used on the same command line with **-rm**
and **-rs** commands.

**Warping meshes and point sets (-rs)**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**-rs <input\_mesh> <output\_mesh>**

The transform chain specified with -r can be applied to points in a
mesh. However, whereas image intensities are mapped from moving space
into fixed space, coordinates of vertices are mapped from fixed space to
moving space.

> greedy -d 3 \\

-rf fixed.nii.gz \\

-rm moving.nii.gz resliced.nii.gz \\

-rs fixed\_mesh.vtk output\_mesh.vtk \\

-r warp1.nii.gz warp2.nii.gz affine.mat

**Matching by Moments of Inertia**
----------------------------------

**Moments mode (-moments)**
~~~~~~~~~~~~~~~~~~~~~~~~~~~

**-moments <1\|2>**

**-o <moments\_matrix>**

Matching by moments can be an effective strategy when two images have
similar content, but are so misaligned that the affine and rigid modes
fail. Matching by moments is particularly useful for binary objects,
e.g. two hippocampal segmentations. Matching by moments line up the
centers of mass of the two images, and (optionally) match the second
momentum tensors.

-  If the argument to **-moments** is 1, only centers of mass are
       matched.

-  If the argument to **-moments** is 2, the second moment tensors are
       also matched.

   -  These is ambiguity with respect to reflection when matching the
          second tensors. Greedy will consider all possible reflections
          (e.g., in 3D there are 8 possible reflections) and choose the
          one that minimizes the metric between fixed and moving images.

The output is a matrix file, just as in affine and rigid registration.
However, unlike rigid and affine modes, the matrix may also include a
coordinate flip (reflection).

**Restricting flipping in moments mode (-det)**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**-det <1\|-1>**

For a 3D image there are 8 ways to line up second momentum tensors,
since the direction is along each momentum vector is arbitrary. Four of
these ways involve flipping, and four do not. By default, the alignment
of tensors that gives the best metric value is used. However, you can
force flipping to always occur (e.g., when you know that you are
matching a left hippocampus mask to a right hippocampus mask) by setting
**-det -1**. Likewise you can prevent flipping by setting **-det 1**.
This option has no effect when using **-moments 1**.

**Disabling rotation in moments mode (-cov-id)**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**-cov-id**

This option sets the second moment tensors to have identity covariance,
which means the matching will not perform any rotation, only alignment
of centers of mass and flipping. Note that **-moments 2 -cov-id** will
allow flipping, whereas **-moments 1** only aligns the centers of mass.
