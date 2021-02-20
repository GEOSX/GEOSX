.. _Continuous_Integration_process:

Continuous Integration process
==============================

To save building time, the third party libraries (that do not change so often) and GEOSX are build separately.

Everytime a pull is requested in the TPL repository, docker images are generated and deployed on `dockerhub <https://hub.docker.com/r/geosx>`_.
The repository names (`ubuntu18.04-gcc8 <https://hub.docker.com/r/geosx/ubuntu18.04-gcc8>`_,
`centos7.7.1908-clang9.0.0 <https://hub.docker.com/r/geosx/centos7.5.1804-clang6.0.1>`_, `centos7.6.1810-gcc8.3.1-cuda10.1.243 <https://hub.docker.com/r/geosx/centos7.6.1810-gcc8.3.1-cuda10.1.243>`_ etc.)
obviously reflect the OS and the compiler flavour used.
For each image, the unique ``${TRAVIS_PULL_REQUEST}-${TRAVIS_BUILD_NUMBER}`` tag is used so we can connect the related code source in a rather convenient way.
Each docker contains the ``org.opencontainers.image.created`` and ``org.opencontainers.image.revision`` labels to provide additional information.

For the OSX builds, we construct a tarball of the TPLs and save them in a remote cloud storage.
There is currently only one mac osx tested environment (xcode 11.2) and the same ``${TRAVIS_PULL_REQUEST}-${TRAVIS_BUILD_NUMBER}`` pattern is used as an identifier for the build.
An important counterpart to using a tarball and not a docker image is that the tarball does not provide the whole system the precompiled binaries rely on.
Problems may arise since we use the rolling release `Homebrew <https://brew.sh/>`_ (to install open-mpi in particular).
To circumvent this potential issue, the brew version is fixed to a specific commit (see BREW_HASH variable in `third party's .travis.yml <https://github.com/GEOSX/thirdPartyLibs/blob/master/.travis.yml>`_)
and stored as a metainformation of the tarball blob inside the cloud storage.
It is therefore possible for GEOSX to recover this informatiom and build against the same revision of brew packages.
Note that the ``TRAVIS_PULL_REQUEST``, ``TRAVIS_BUILD_NUMBER`` and ``TRAVIS_COMMIT`` are also stored as metainformation in the same way
(have a look at the OSX build section of `GEOSX's .travis.yml <https://github.com/GEOSX/GEOSX/blob/develop/.travis.yml>`_ to see how to retrieve these informations).

There thus is only one unique identifier for both dockers and mac osx builds for one TPL code base.
It is necessary to define the global environment ``GEOSX_TPL_TAG`` (`e.g.` something like ``82-254``) to build against one selected version of the TPL.

It must be mentioned that one and only one version of the compiled TPL tarball is stored per pull request (older ones are removed automatically).
Therefore, a client building against a work in progress PR may experience a 404 error sooner or later.

Building docker images
----------------------

Our continuous integration process builds the TPL and GEOSX against two operating systems (ubuntu and centos) and two compilers (clang and gcc).
The docker files use `multi-stage builds <https://docs.docker.com/develop/develop-images/multistage-build/>`_ in order to minimise the sizes of the images.

* First stage installs and defines all the elements that are commons to both TPL and GEOSX (for example, MPI and c++ compiler, BLAS, LAPACK, path to the installation directory...).
* As a second stage, we install everything needed to build (`not run`) the TPLs.
  We keep nothing from this second step for GEOSX, except the compiled TPL themselves.
  For example, a fortran compiler is needed by the TPL but not by GEOSX: it shall be installed during this step, so GEOSX won't access a fortran compiler (it does not have to).
* Last stage copies the compiled TPL from second stage and installs the elements only required by GEOSX (there are few).

.. _Docker_images_contract:

Docker images contract
----------------------

GEOSX will find a compiled version of the third party libraries.

As part of the contract provided by the TPL, the docker images also defines several environment variables.
The

.. code-block:: sh

    GEOSX_TPL_DIR

variable contains the absolute path of the installation root directory of the third party libraries.
GEOSX must use it when building.

Other variables are classical absolute path compiler variables.

.. code-block:: sh

    CC
    CXX
    MPICC
    MPICXX

And the absolute path the mpirun (or equivalent) command.

.. code-block:: sh

    MPIEXEC

The following ``openmpi`` environment variables allow it to work properly in the docker container.
But there should be no reason to access or use them explicitly.

.. code-block:: sh

    OMPI_CC=$CC
    OMPI_CXX=$CXX
