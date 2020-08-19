************
Installation
************

Compiling from Source Code
**************************

Building from source code is the best way to obtain the latest version of Greedy and it works on all platforms. The instructions below are for Linux and MacOS.

Prerequisites
=============

1. `ITK`_ 4.12.2 or later 
2. `CMake`_ 3.9 or later (and basic familiarity with this tool)
3. `Git`_ (and basic familiarity with this tool)

Building Greedy
===============

Set up the directory tree and clone greedy repository::

    cd my_directory
    git clone https://github.com/pyushkevich/greedy greedy
    mkdir build

At this point, ``my_directory/greedy`` contains the source code and ``my_directory/build`` is where we will compile greedy.

Run ``ccmake`` from the build directory::

    cd build
    ccmake ../greedy

When running ``ccmake``, set the following variables:

==================    ===================
Variable              Value
==================    ===================
ITK_DIR               Point to the directory where you compiled ITK
CMAKE_BUILD_TYPE      Release (unless you want to debug)
USE_FFTW              OFF (you can set it to on, but it's not going to affect greedy)
==================    ===================

After you configure and generate makefiles in ``ccmake``, run::

    make

You should now find the executable ``greedy`` created in the ``build`` directory. You can test that Greedy was compiled by running::

    ./greedy -v


Using Pre-compiled Binaries
***************************

As of late 2018, greedy is part of `ITK-SNAP`_ software. Simply download ITK-SNAP 3.8 or later for your platform and run ``Help->Install Command-Line Tools`` from the ITK-SNAP main menu.


.. _ITK: http://itk.org/
.. _CMake: http://cmake.org/
.. _Git: https://git-scm.com/
.. _ITK-SNAP: http://itksnap.org
