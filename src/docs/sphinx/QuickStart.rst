###############################
Quick Start Guide
###############################

The goal of this page is to get you started as quickly as possible using GEOSX.  We will walk you through downloading the source, compiling the code, and testing the installation.

Before jumping to the installation process, we want to first address some frequently asked questions we get from new users.  If you are itching to get started, feel free to jump ahead to the relevant sections.

Frequently Asked Questions
==========================

Does GEOSX have a graphical user interface?
------------------------------------------------
Given the focus on rapid development and HPC environments, GEOSX does not have a graphical user interface.
This is consistent with many other high performance computing packages, but we recognize it can be a deal-breaker for certain users.
For those who can get past this failing, we promise we still have a lot to offer.
In a typical workflow, you will prepare an XML-based input file describing your problem.
You may also prepare a mesh file containing geometric and property information describing, say, a reservoir you would like to simulate.
There is no shortage of GUI tools that can help you in this model building stage.
The resulting input deck is then consumed by GEOSX to run the simulation and produce results.
This may be done in a terminal of your local machine or by submitting a job to a remote server.
The resulting output files can then be visualized by any number of graphical visualization programs (typically `VisIt <https://wci.llnl.gov/simulation/computer-codes/visit/>`_ or `paraview <https://www.paraview.org/>`_).
Thus, while GEOSX is GUI free, the typical workflow is not.

Do I need to be a code developer to use GEOSX?
------------------------------------------------
For the moment, most users will
need to download and compile the code from source, which we readily admit this requires
a certain level of development expertise.  We try to make this process as easy as
possible, and we are working on additional deployment options to make this process easier.  Once installed, however, our goal is to make GEOSX accessible to developers and non-developers alike.  Our target audience includes engineers and scientists who want to solve tough application problems, but could care less about the insides of the tool.  For those of you who *are* interested in scientific computing, however, GEOSX is an open source project and we welcome external contributions.

What are the system requirements?
------------------------------------------------
GEOSX is primarily written in C++, with a focus on standards compliance and platform-to-platform portability.
It is designed to run on everything from commodity laptops to the world's most powerful supercomputers.
We regularly test the code across a variety of operating systems and compilers.
Most of these operating systems are Linux/UNIX based (e.g. Ubuntu, CentOS, Mac OSX).
We do have developers working in Windows environments, but they use a Virtual Machine or work within a docker image rather than directly in the Windows environment.
In the instructions below, we assume you have access to fairly standard development tools.
Using advanced features of GEOSX, like GPU-acceleration, will of course introduce additional hardware and software requirements.

Help, I get errors while trying to download/compile/run!
---------------------------------------------------------
Unfortunately, no set of instructions is foolproof.  It is simply impossible to anticipate every system configuration or user.  If you run into problems during the installation, we recommend the following five-step process:

#. Take a moment to relax, and then re-read the instructions carefully.  Perhaps you overlooked a key step?  Re-read the error message(s) closely.  Modern compilation tools are often quite helpful in reporting exactly why things fail.

#. Type a few keywords from your error into a search engine.  It is possible someone else out there has encountered your problem before, and a well-chosen keyword can often produce an instant solution.  Note that when a compilation fails, you may get pages and pages of errors.  Try to identify the *first* one to occur and fix that.  One error will often trigger subsequent errors, and looking at the *last* error on the screen may not be so helpful.

#. If you encounter problems building one of the third-party libraries we depend on, check out their support pages.  They may be able to help you more directly than we can.

#. Still stuck? Check out our `issues tracker <https://github.com/GEOSX/GEOSX/issues>`_, searching current or closed issues that may address your problem.  Perhaps someone has had an identical issue, or something close.  The issue tracker has a convenient search bar where you can search for relevant keywords.  Remember to remove the default ``is:open`` keyword to search both open and closed issues.

#. If you have exhausted the options above, it is time to seek help from the developers.  Post an issue on our issue tracker.  Be specific, providing as much information as possible about your system setup and the error you are encountering.  Please be patient in this process, as we may need to correspond a few times and ask you to run additional tests.  Most of the time, users have a slightly unusual system configuration that we haven't encountered yet, such as an older version of a particular library.  Other times there is a legitimate bug in GEOSX to be addressed.  Take pride in the fact that you may be saving the next user from wasted time and frustration.

Repository Organization
==============================
The source for GEOSX and related tools are hosted on `Github <https://github.com>`_.
We use `Git workflows <https://git-scm.com>`_ to version control our code and manage the entire development process.
On Github, we have a `GEOSX Organization <https://github.com/GEOSX>`_ that hosts several related repositories.

You should sign up for a free Github account, particularly if you are interested in posting issues to our issue tracker and communicating with the developers. The main repository of interest is obviously the `GEOSX repository <https://github.com/GEOSX/GEOSX>`_ itself.

We also rely on two types of dependencies: first-party and third-party.
First-party dependencies are projects directly associated with the GEOSX effort, but kept in separate repositories because they form stand-alone tools.
For example, there is a geologic mesh handling package called `PAMELA <https://github.com/GEOSX/PAMELA>`_ and an equation-of-state package called `PVTPackage <https://github.com/GEOSX/PVTPackage>`_.
These packages are handled as `Git Submodules <https://git-scm.com/book/en/v2/Git-Tools-Submodules>`_, which provides a transparent way of coordinating multiple code development projects.
Most users will never have to worry that these modules are in fact separate projects from GEOSX.

We also rely on several open-source Third-Party Libraries (TPLs) (see `thirdPartyLibs <https://github.com/GEOSX/thirdPartyLibs>`_).
These are well-respected projects developed externally to GEOSX.  We have found, however, that many compilation issues stem from version incompatibilities between different packages.
To address this, we provide a mirror of these TPLs, with version combinations we know play nicely together.
We also provide a build script that conveniently and consistently builds those dependencies.

Our build system will automatically use the mirror package versions by default.
You are welcome to tune your configuration, however, to point to different versions installed on your system.
If you work on an HPC platform, for example, common packages may already be available and optimized for platform hardware.
For new users, however, it may be safer to begin with the TPL mirror.

.. note::
   If you are working on an HPC platform with several other GEOSX users, we often compile the TPLs in a shared location so individual users don't have to waste their storage quota. Inquire with your institution's point-of-contact whether this option already exists. For all LLNL systems, the answer is yes.

Finally, there are also several private repositories only accessible to the core development team, which we use for behind-the-scene testing and maintenance of the code.

Username and Authentication
=============================
New users should sign up for a free `Github account <https://github.com>`_.
If you intend to develop in the GEOSX codebase, you may benefit from setting up your git credentials (see :ref:`GitWorkflow`).


Download
======================

It is possible to directly download the source code as a zip file.
We strongly suggest, however, that users don't rely on this option.
Instead, most users should use Git to either *clone* or *fork* the repository.
This makes it much easier to stay up to date with the latest releases and bug fixes.
If you are not familiar with the basics of Git, `here is a helpful resource <https://git-scm.com>`_ to get you started.

The tutorial here assumes you will use a https clone with no specific credentials.
Using an ssh connection pattern requires a very slight modification (described in the notes below).

If you do not already have Git installed on your system, you will need to install it.
We recommend using a relatively recent version of Git, as there have been some notable improvements over the past few years.
You can check if Git is already available by opening a terminal and typing

.. code-block:: sh

  git --version

You'll also need the `git-lfs <https://git-lfs.github.com/>`_ large file extension.

The first task is to clone the ``GEOSX`` and ``thirdPartyLibs`` repositories.
If you do not tell it otherwise, the build system will expect the GEOSX and thirdPartyLibs to be parallel to each other in the directory structure.  For example,

.. code-block:: sh

  codes/
  ├── GEOSX/
  └── thirdPartyLibs/

where the toplevel ``codes`` directory can be re-named and located wherever you like. It is possible to customize the build system to expect a different structure, but for now let us assume you take the simplest approach.

First, using a terminal, create the ``codes`` directory wherever you like.

.. code-block:: sh

  cd /insert/your/desired/path/
  mkdir codes
  cd codes

Inside this directory, we can clone the GEOSX repository.  We will also use some Git commands to initialize and download the submodules (e.g. PAMELA). Note that most users will not have access to our integrated tests repository, and so we
"deinit" (deactivate) this submodule.  Developers who will be working with the integratedTests repository should skip this line.

.. code-block:: sh

   git clone https://github.com/GEOSX/GEOSX.git
   cd GEOSX
   git lfs install
   git submodule init
   git submodule deinit integratedTests
   git submodule update
   cd ..

If all goes well, you should have a complete copy of the GEOSX source at this point.
The most common errors people encounter here have to do with Github not recognizing their authentication settings and/or repository permissions.


.. note::
   To use SSH instead, the clone step becomes ``git clone git@github.com:GEOSX/GEOSX.git``. See the previous section for tips on ensuring your SSH is working properly.  Alternatively, you may want to insert credentials in the command line (though less secure) ``git clone https://${USER}:${TOKEN}@github.com/GEOSX/GEOSX.git``.

.. note::
   The integratedTests submodule is not publicly available, with access limited to the core development team.  This may cause the ``git submodule update`` command to fail
   if you forget the ``git submodule deinit integratedTests`` step above.
   This submodule is not required for building GEOSX. If you see an error message here, however, you may need to initialize and update the submodules manually:

   .. code-block:: sh

      cd GEOSX
      git submodule init src/coreComponents/LvArray
      git submodule init src/coreComponents/fileIO/coupling/hdf5_interface
      git submodule init src/coreComponents/physicsSolvers/GEOSX_PTP
      git submodule init src/externalComponents/PAMELA
      git submodule init src/externalComponents/PVTPackage
      git submodule update src/coreComponents/LvArray
      git submodule update src/coreComponents/fileIO/coupling/hdf5_interface
      git submodule update src/coreComponents/physicsSolvers/GEOSX_PTP
      git submodule update src/externalComponents/PAMELA
      git submodule update src/externalComponents/PVTPackage
      cd ..

.. note::
   ``git-lfs`` may not function properly (or may be very slow) if your version of git and git-lfs are not current. If you are using an older version, you may need to add ``git lfs pull`` after ``git pull`` in the above procedures.

Once we have grabbed GEOSX, we do the same for the thirdPartyLibs repository.  From the ``codes`` directory, type

.. code-block:: sh

   git clone https://github.com/GEOSX/thirdPartyLibs.git
   cd thirdPartyLibs
   git lfs install
   git pull
   git submodule init
   git submodule update
   cd ..

Again, if all goes well you should now have a copy of all necessary TPL packages.


Configuration
================

At a minimum, you will need a relatively recent compiler suite installed on your system (e.g. `GCC <https://gcc.gnu.org>`_, `Clang <https://clang.llvm.org>`_) as well as `CMake <https://cmake.org>`_.  If you want to run jobs using MPI-based parallelism, you will also need an MPI implementation (e.g. `OpenMPI <https://www.open-mpi.org>`_, `MVAPICH <https://mvapich.cse.ohio-state.edu>`_).  Note that GEOSX supports a variety of parallel computing models, depending on the hardware and software environment.  Advanced users are referred to the :ref:`BuildGuide` for a discussion of the available configuration options.

Before beginning, it is a good idea to have a clear idea of the flavor and version of the build tools you are using.  If something goes wrong, the first thing the support team will ask you for is this information.

.. code-block:: sh

  cpp --version
  mpic++ --version
  cmake --version

Here, you may need to replace ``cpp`` with the full path to the C++ compiler you would like to use, depending on how your path and any aliases are configured.

GEOSX compilations are driven by a cmake ``host-config`` file, which tells the build system about the compilers you are using, where various packages reside, and what options you want to enable.  We have created a number of default hostconfig files for common systems.  You should browse them to see if any are close to your needs:

.. code-block:: sh

   cd GEOSX/host-configs

We maintain host configs (ending in ``.cmake``) for HPC systems at various institutions, as well as ones for common personal systems.  If you cannot find one that matches your needs, we suggest beginning with one of the shorter ones and modifying as needed.  A typical one may look like:

.. code-block:: sh

  # file: your-platform.cmake

  # detect host and name the configuration file
  site_name(HOST_NAME)
  set(CONFIG_NAME "your-platform" CACHE PATH "")
  message( "CONFIG_NAME = ${CONFIG_NAME}" )

  # set paths to C, C++, and Fortran compilers. Note that while GEOSX does not contain any Fortran code,
  # some of the third-party libraries do contain Fortran code. Thus a Fortran compiler must be specified.
  set(CMAKE_C_COMPILER "/usr/bin/clang" CACHE PATH "")
  set(CMAKE_CXX_COMPILER "/usr/bin/clang++" CACHE PATH "")
  set(CMAKE_Fortran_COMPILER "/usr/local/bin/gfortran" CACHE PATH "")
  set(ENABLE_FORTRAN OFF CACHE BOOL "" FORCE)

  # enable MPI and set paths to compilers and executable.
  # Note that the MPI compilers are wrappers around standard serial compilers.
  # Therefore, the MPI compilers must wrap the appropriate serial compilers specified
  # in CMAKE_C_COMPILER, CMAKE_CXX_COMPILER, and CMAKE_Fortran_COMPILER.
  set(ENABLE_MPI ON CACHE PATH "")
  set(MPI_C_COMPILER "/usr/local/bin/mpicc" CACHE PATH "")
  set(MPI_CXX_COMPILER "/usr/local/bin/mpicxx" CACHE PATH "")
  set(MPI_Fortran_COMPILER "/usr/local/bin/mpifort" CACHE PATH "")
  set(MPIEXEC "/usr/local/bin/mpirun" CACHE PATH "")

  # disable CUDA and OpenMP
  set(CUDA_ENABLED "OFF" CACHE PATH "" FORCE)
  set(ENABLE_OPENMP "OFF" CACHE PATH "" FORCE)

  # enable PAMELA and PVTPackage
  set(ENABLE_PAMELA ON CACHE BOOL "" FORCE)
  set(ENABLE_PVTPackage ON CACHE BOOL "" FORCE)

  # enable tests
  set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "" FORCE )

The various ``set()`` commands are used to set environment variables that control the build.  You will see in the above example that we set the C++ compiler to ``/user/bin/clang++`` and so forth.  We also disable CUDA and OpenMP, but enable PAMELA and PVTPackage.  The final line is related to our unit test suite.  See the :ref:`BuildGuide` for more details on available options.

.. note::
   If you develop a new ``host-config`` for a particular platform that may be useful for other users, please consider sharing it with the developer team.

Compilation
==================

We will begin by compiling the TPLs, followed by the main code.  If you work on an HPC system with other GEOSX developers, check with them to see if the TPLs have already been compiled in a shared directory.  If this is the case, you can skip ahead to just compiling the main code.  If you are working on your own machine, you will need to compile both.

We strongly suggest that GEOSX and TPLs be built with the same hostconfig file.  Below, we assume that you keep it in, say, ``GEOSX/host-configs/your-platform.cmake``, but this is up to you.

We begin with the third-party libraries, and use a python ``config-build.py`` script to configure and build all of the TPLs.  Note that we will request a Release build type, which will enable various optimizations.  The other option is a Debug build, which allows for debugging but will be much slower in production mode.  The TPLS will then be built in a build directory named consistently with your hostconfig file.

.. code-block:: sh

   cd thirdPartyLibs
   python scripts/config-build.py -hc ../GEOSX/host-configs/your-platform.cmake -bt Release
   cd build-your-platform-release
   make

Note that building all of the TPLs can take quite a while, so you may want to go get a cup of coffee at this point.
Also note that you should *not* use a parallel ``make -j N`` command to try and speed up the build time.

The next step is to compile the main code.
Again, the ``config-build.py`` sets up cmake for you, so the process is very similar.

.. code-block:: sh

   cd ../../GEOSX
   python scripts/config-build.py -hc host-configs/your-platform.cmake -bt Release
   cd build-your-platform-release
   make -j4
   make install

Here, the parallel ``make -j 4`` will use four processes for compilation, which can substantially speed up the build if you have a multi-processor machine.
You can adjust this value to match the number of processors available on your machine.
The ``make install`` command then installs GEOSX to a default location unless otherwise specified.
If all goes well, a ``geosx`` executable should now be available:

.. code-block:: sh

  GEOSX/install-your-platform-release/bin/geosx

Running
=================

We can do a quick check that the geosx executable is working properly by calling the executable with our help flag

.. code-block:: sh

  ./bin/geosx --help

This should print out a brief summary of the available command line arguments:

.. code-block:: sh

    USAGE: geosx -i input.xml [options]

    Options:
    -?, --help
    -i, --input,            Input xml filename (required)
    -r, --restart,          Target restart filename
    -x, --x-partitions,     Number of partitions in the x-direction
    -y, --y-partitions,     Number of partitions in the y-direction
    -z, --z-partitions,     Number of partitions in the z-direction
    -s, --schema,           Name of the output schema
    -b, --use-nonblocking,  Use non-blocking MPI communication
    -n, --name,             Name of the problem, used for output
    -s, --suppress-pinned   Suppress usage of pinned memory for MPI communication buffers
    -o, --output,           Directory to put the output files
    -t, --timers,           String specifying the type of timer output.
    An input xml must be specified!

Obviously this doesn't do much interesting, but it will at least confirm that the executable runs.  In typical usage, an input XML must be provided describing the problem to be run, e.g.

.. code-block:: sh

    ./bin/geosx -i your-problem.xml

In a parallel setting, the command might look something like

.. code-block:: sh

    mpirun -np 8 ./bin/geosx -i your-problem.xml -x 2 -y 2 -z 2

Note that we provide a series of :ref:`Tutorials` to walk you through the actual usage of the code, with several input examples.  Once you are comfortable the build is working properly, we suggest new users start working through these tutorials.

Testing
=================

It is wise to run our unit test suite as an additional check that everything is working properly. You can run them in the build folder you just created.

.. code-block:: sh

  cd GEOSX/build-your-platform-release
  ctest -V

This will run a large suite of simple tests that check various components of the code.  If you have access, you may also consider running the integrated tests.  Please refer to :ref:`IntegratedTests` for further information.

.. note::
   If *all* of the unit tests fail, there is likely something wrong with your installation.  Refer to the FAQs above for how best to proceed in this situation.  If only a few tests fail, it is possible that your platform configuration has exposed some issue that our existing platform tests do not catch.  If you suspect this is the case, please consider posting an issue to our issue tracker (after first checking whether other users have encountered a similar issue).
