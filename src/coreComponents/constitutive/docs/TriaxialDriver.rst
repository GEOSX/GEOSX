Triaxial Driver
===============

When calibrating solid material parameters to experimental data, it can be a hassle to launch a full finite element simulation to mimic experimental loading conditions.  Instead, GEOSX provides a ``TriaxialDriver`` allowing the user to run loading tests on a single material point.  This makes it easy to understand the material response and fit it to lab data.  The driver itself is launched like any other GEOSX simulation, but with a particular XML structure:

.. code-block:: sh

   ./bin/geosx -i myTest.xml

XML Structure
-------------
A typical XML file to run the triaxial driver will have the following key elements.  We present the whole file first, before digging into the individual blocks.

.. code-block:: xml

    <Problem>

      <!-- Triaxial driver is added as an executable Task-->
      <Tasks>
        <TriaxialDriver
          name="triaxialDriver"
          material="sand"
          mode="triaxial"
          strainFunction="strainFunction"
          stressFunction="stressFunction"
          steps="40"
          output="results.txt"
          logLevel="1" />
      </Tasks>

      <!-- This Task is added to the Event queue as a SoloEvent-->
      <Events
        maxTime="1">
        <SoloEvent
          name="triaxialDriver"
          target="/Tasks/triaxialDriver"/>
      </Events>

      <!-- The driver calls the material "sand" which is defined here-->
      <Constitutive>
        <ExtendedDruckerPrager
          name="sand"
          defaultDensity="2700"
          defaultBulkModulus="500"
          defaultShearModulus="300"
          defaultCohesion="0.0"
          defaultInitialFrictionAngle="15"
          defaultResidualFrictionAngle="23"
          defaultDilationRatio="1.0"
          defaultHardening="0.001"
        />
      </Constitutive>

      <!-- The axial/radial loading conditions are defined by time-dependent functions-->
      <Functions>
        <TableFunction
          name="strainFunction"
          inputVarNames="{ time }"
          coordinates="{ 0.0, 3.0, 4.0, 7.0, 8.0 }"
          values="{ 0, -0.003, -0.002, -0.005, -0.004 }"/>
        <TableFunction
          name="stressFunction"
          inputVarNames="{ time }"
          coordinates="{ 0.0, 8.0  }"
          values="{ -1.0, -1.0 }"/>
      </Functions>

      <!-- A mesh is not actually used, but GEOSX throws an error without one.  Will fix this soon-->
      <Mesh>
        <InternalMesh
          name="mesh1"
          elementTypes="{ C3D8 }"
          xCoords="{ 0, 1 }"
          yCoords="{ 0, 1 }"
          zCoords="{ 0, 1 }"
          nx="{ 1 }"
          ny="{ 1 }"
          nz="{ 1 }"
          cellBlockNames="{ cellBlock01 }"/>
      </Mesh>

    </Problem>

The first thing to note is that the XML structure is identical to a standard GEOSX input deck.  In fact, once the constitutive block is calibrated, one could start adding solver and discretization blocks to the same file to create a proper field simulation.  This makes it easy to go back and forth between calibration and simulation.

The ``TriaxialDriver`` is added as a ``Task``, a particular type of executable event often used for simple actions.  It is added as a ``SoloEvent`` to the event queue.  This leads to a trivial event queue, since all we do is launch the driver and then quit.

.. note::

   Internally, the triaxial driver uses a simple form of time-stepping to advance through the loading steps, allowing for both rate-dependent and rate-independent models to be tested. This timestepping is handled independently from the more complicated time-stepping pattern used by physics ``Solvers`` and coordinated by the ``EventManager``.  In particular, in the XML file above, the ``maxTime`` parameter in the ``Events`` block is an event manager control, controlling when/if certain events occur.  Once launched, the triaxial driver internally determines its own max time and timestep size using a combination of the strain function's time coordinates and the requested number of loadsteps.  It is therefore helpful to think of the driver as an instantaneous *event* (from the event manager's point of view), but one which has a separate, internal clock.

The key parameters for the TriaxialDriver are:

.. include:: /coreComponents/fileIO/schema/docs/TriaxialDriver.rst

.. note::

   GEOSX uses the *engineering* sign convention where compressive stresses and strains are *negative*.
   This is one of the most frequent issues users make when calibrating material parameters, as
   stress- and strain-like quantities often need to be negative to make physical sense.  You may note in the
   XML above, for example, that ``stressFunction`` and ``strainFunction`` have negative values for
   a compressive test.

Test Modes
----------
The most complicated part of the driver is understanding how the stress and strain functions are applied in different testing modes.  The driver mimics laboratory core tests, with loading controlled in the
axial and radial directions. These conditions may be either strain-controlled or stress-controlled.  The following table describes the available test modes in detail:

+---------------+-------------------------+-------------------------+---------------------------+
| **mode**      | **axial loading**       | **radial loading**      | **initial stress**        |
+---------------+-------------------------+-------------------------+---------------------------+
| ``triaxial``  | axial strain controlled | radial stress controlled| isotropic stress using    |
|               | with ``strainFunction`` | with ``stressFunction`` | ``stressFunction(t=tmin)``|
+---------------+-------------------------+-------------------------+---------------------------+
| ``volumetric``| axial strain controlled | radial strain =         | isotropic stress using    |
|               | with ``strainFunction`` | axial strain            | ``stressFunction(t=tmin)``|
+---------------+-------------------------+-------------------------+---------------------------+
| ``oedometer`` | axial strain controlled | zero radial strain      | isotropic stress using    |
|               | with ``strainFunction`` |                         | ``stressFunction(t=tmin)``|
+---------------+-------------------------+-------------------------+---------------------------+

To set the initial stress state, the ``stressFunction`` is evaluated at ``t=tmin`` (usually ``t=0``, though conceivably a user may put in a time-function with a non-zero starting point). This scalar value is used to set the material to an isotropic initial stress state.  In the volumetric and oedometer tests, the remainder of the ``stressFunction`` time history is ignored, as they are strain-controlled tests.  By setting the initial stress this way, it makes it easy to start the test from a well-defined confining pressure.

Output Format
-------------
The ``output`` key is used to identify a file to which the results of the simulation are written.  If this key is omitted, or the user specifies ``output="none"``, file output will be suppressed.  The file is a simple ASCII format with a brief header followed by test data:

.. code:: sh

  # column 1 = time
  # column 2 = axial_strain
  # column 3 = radial_strain_1
  # column 4 = radial_strain_2
  # column 5 = axial_stress
  # column 6 = radial_stress_1
  # column 7 = radial_stress_2
  # column 8 = newton_iter
  # column 9 = residual_norm
  0.0000e+00  0.0000e+00 0.0000e+00 0.0000e+00 -1.0000e+00 -1.0000e+00 -1.0000e+00 0.0000e+00 0.0000e+00
  1.6000e-01 -1.6000e-04 4.0000e-05 4.0000e-05 -1.1200e+00 -1.0000e+00 -1.0000e+00 2.0000e+00 0.0000e+00
  3.2000e-01 -3.2000e-04 8.0000e-05 8.0000e-05 -1.2400e+00 -1.0000e+00 -1.0000e+00 2.0000e+00 0.0000e+00
  ...

This file can be readily plotted using any number of plotting tools.  Each row corresponds to one timestep of the driver, starting from initial conditions in the first row.

We note that the file contains two columns for radial strain and two columns for radial stress.  For an isotropic material, the stresses and strains along the two radial axes will usually be identical.  We choose to output this way, however, to accommodate both anisotropic materials and true-triaxial loading conditions.  In these cases, the stresses and strains in the radial directions could potentially differ.

These columns can be added and subtracted to produce other quantities of interest, like mean stress or deviatoric stress.  For example, we can plot the output of our "sand" XML to produce the following stress / strain curves:

.. figure:: TriaxialDriver.svg
   :width: 600px
   :align: center
   :alt: stress/strain figure
   :figclass: align-center

   **Figure**: Stress/strain behavior resulting from the sand model XML above.

In this plot, we have reversed the sign convention to be consistent with typical experimental plots.  Note also that the ``strainFunction`` above includes two unloading cycles, allowing us to observe both plastic loading and elastic unloading.

Model Convergence
-----------------

The last two columns of the output file contain information about the convergence behavior of the material driver.  In ``triaxial`` mode, the mixed nature of the stress/strain control requires using a Newton solver to converge the solution.  This last column reports the number of Newton iterations and final residual norm.  Large values here would be indicative of the material model struggling (or failing) to converge.  Convergence failures can result from several reasons, including:

1. Inappropriate material parameter settings
2. Overly large timesteps
3. Infeasible loading conditions (i.e. trying to load a material to a physically-unreachable stress point)
4. Poor model implementation

We generally spend a lot of time vetting the material model implementations (#4).  When you first encounter a problem, it is therefore good to explore the other three scenarios first.  If you find something unusual in the model implementation or are just really stuck, please submit an issue on our issue tracker so we can help resolve any bugs.

Unit Testing
------------

The development team also uses the Triaxial Driver to perform unit testing on the various material models within GEOSX.  The optional argument ``baseline`` can be used to point to a previous output file that has been validated  (e.g. against analytical or experimental benchmarks).  If such a file is specified, the driver will perform a loading run and then compare the new results against the baseline.  In this way, any regressions in the material models can be quickly identified.  

Developers of new models are encouraged to add their own baselines to ``src/coreComponents/constitutive/unitTests``. Adding additional tests is straightforward:

1. Create a new xml file for your test in ``src/coreComponents/constitutive/unitTests``.  There are several examples is this directory already to use as a template.  We suggest using the naming convention ``testTriaxial_myTest.xml``, so that all triaxial tests will be grouped together alphabetically.  Set the ``output`` file to ``testTriaxial_myTest.txt``, and run your test.  Validate the results however is appropriate.  

2. This output file will now become your new baseline.  Replace the ``output`` key with ``baseline`` so that the driver can read in your file as a baseline for comparison.  Make sure there is no remaining ``output`` key, or set ``output=none``, to suppress further file output.  While you can certainly write a new output for debugging purposes, during our automated unit tests we prefer to suppress file output.  Re-run the triaxial driver to confirm that the comparison test passes.

3. Modify ``src/coreComponents/constitutive/unitTests/CMakeLists.txt`` to enable your new test in the unit test suite.  In particular, you will need to add your new XML file to the existing list in the ``gtest_triaxial_xmls`` variable:

.. code:: sh

  set( gtest_triaxial_xmls
       testTriaxial_elasticIsotropic.xml
       testTriaxial_druckerPragerExtended.xml
       testTriaxial_myTest.xml
     )
4. Run ``make`` in your build directory to make sure the CMake syntax is correct

5. Run ``ctest -V -R Triax`` to run the triaxial unit tests.  Confirm your test is included and passes properly.

If you run into troubles, do not hesitate to contact the development team for help.



