/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ProblemManager.hpp
 */


#ifndef GEOSX_MANAGERS_PROBLEMMANAGER_HPP_
#define GEOSX_MANAGERS_PROBLEMMANAGER_HPP_

#ifdef GEOSX_USE_PYTHON
// Note: the python header must be included first to avoid conflicting
// definitions of _posix_c_source
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>
#endif

#include "EventManager.hpp"
#include "managers/Functions/FunctionManager.hpp"
#include "mesh/GraphManager.hpp"
#include "fileIO/schema/schemaUtilities.hpp"

namespace geosx
{

class PhysicsSolverManager;
class DomainPartition;
namespace constitutive
{
class ConstitutiveManager;
}
/**
 * @class ProblemManager
 * @brief This is the class handling the operation flow of the problem being ran in GEOSX
 */
class ProblemManager : public dataRepository::Group
{
public:

  /**
   * @param name the name of this object manager
   * @param parent the parent Group
   */
  explicit ProblemManager( const std::string & name,
                           Group * const parent );

  /**
   * @brief Destructor, deletes all Groups and Wrappers owned by this Group
   */
  ~ProblemManager() override;

  /**
   * @brief Handles deviations between the data structure and schema
   * @param schemaRoot schema root node handle
   * @param schemaParent schema parent node handle
   * @param documentationType flag to indicate the type of schema (0=input, 1=other)
   * @details This function handles deviations between the xml and data structure
   * on the Problem level (Functions, Mesh, etc.).  This can also be used to
   * add entries to the schema, which are not used during normal code execution
   * (e.g.: Benchmark)
   */
  virtual void SetSchemaDeviations( xmlWrapper::xmlNode schemaRoot,
                                    xmlWrapper::xmlNode schemaParent,
                                    integer documentationType ) override;

  /**
   * @brief Creates a new sub-Group using the ObjectCatalog functionality.
   * @param childKey The name of the new object type's key in the
   *                 ObjectCatalog.
   * @param childName The name of the new object in the collection of
   *                  sub-Groups.
   * @return A pointer to the new Group created by this function.
   */
  virtual Group * CreateChild( string const & childKey, string const & childName ) override;

  /**
   * @brief Parses command line input
   */
  void ParseCommandLineInput();

  /**
   * @brief Parses a restart file
   * @param restartFileName the name of the restart file
   * @return flag indicating beginFromRestart status
   */
  static bool ParseRestart( std::string & restartFileName );

  /**
   * @brief Initializes a python interpreter within GEOSX
   * @note This is not regularly used or tested, and may be removed in future versions.
   * To use this feature, the code must be compiled with the GEOSX_USE_PYTHON flag
   */
  void InitializePythonInterpreter();

  /**
   * @brief Closes the internal python interpreter
   * @note This is not regularly used or tested, and may be removed in future versions.
   * To use this feature, the code must be compiled with the GEOSX_USE_PYTHON flag
   */
  void ClosePythonInterpreter();

  /**
   * @brief Generates the xml schema documentation
   * This function is called when the code is called with the -s schema_name option.
   * @details Before generating the schema, the code builds up a comprehensive datastructure.
   * (Note: catalog objects throughout the code will typically be registered via the
   * ExpandObjectCatalogs method.)  Once ready, SchemaUtilities will recusively walk
   * through the database, generating the xml schema.
   */
  void GenerateDocumentation();

  /**
   * @brief Parses the input xml file
   * @details The name of the input file is indicated via the -i option on the command line
   */
  void ParseInputFile();

  /**
   * @brief Generates numerical meshes used throughout the code
   */
  void GenerateMesh();

 /**
   * @brief Generates graphs used throughout the code
   */
  void GenerateGraph();


  /**
   * @brief Allocates constitutive relations according to the discretizations
   *   on each subregion.
   */
  void ApplyNumericalMethods();

  /**
   * @brief Defines the order in which objects should be initialized
   * @param order list defining ordering sequence
   */
  void InitializationOrder( string_array & order ) override final;

  /**
   * @brief Sets up the problem after the input has been read in
   */
  void ProblemSetup();

  /**
   * @brief Run the events in the scheduler.
   */
  void RunSimulation();

  /**
   * @brief After initialization, overwrites data using a restart file
   */
  void ReadRestartOverwrite();

  /**
   * @brief Applies initial conditions indicated within the input file FieldSpecifications block
   */
  void ApplyInitialConditions();

  /**
   * @brief Returns a pointer to the DomainPartition
   * @return Pointer to the DomainPartition
   */
  DomainPartition * getDomainPartition();

  /**
   * @brief Returns a pointer to the DomainPartition
   * @return Const pointer to the DomainPartition
   */
  DomainPartition const * getDomainPartition() const;

  /**
   * @brief Returns the problem name
   * @return The problem name
   */
  const string & getProblemName() const
  { return GetGroup< Group >( groupKeys.commandLine )->getReference< string >( viewKeys.problemName ); }

  /**
   * @brief Returns the input file name
   * @return The input file name
   */
  const string & getInputFileName() const
  { return GetGroup< Group >( groupKeys.commandLine )->getReference< string >( viewKeys.inputFileName ); }

  /**
   * @brief Returns the restart file name
   * @return The restart file name
   */
  const string & getRestartFileName() const
  { return GetGroup< Group >( groupKeys.commandLine )->getReference< string >( viewKeys.restartFileName ); }

  /**
   * @brief Returns the schema file name
   * @return The schema file name
   */
  const string & getSchemaFileName() const
  { return GetGroup< Group >( groupKeys.commandLine )->getReference< string >( viewKeys.schemaFileName ); }

  /// Input file xml document handle
  xmlWrapper::xmlDocument xmlDocument;

  /// Input file parsing results
  xmlWrapper::xmlResult xmlResult;

  /// Input file Problem node handle
  xmlWrapper::xmlNode xmlProblemNode;

  /// Command line input viewKeys
  struct viewKeysStruct
  {
    dataRepository::ViewKey inputFileName            = {"inputFileName"};            ///< Input file name key
    dataRepository::ViewKey restartFileName          = {"restartFileName"};          ///< Restart file name key
    dataRepository::ViewKey beginFromRestart         = {"beginFromRestart"};         ///< Flag to begin from restart key
    dataRepository::ViewKey xPartitionsOverride      = {"xPartitionsOverride"};      ///< Override of number of
                                                                                     ///< subdivisions in x key
    dataRepository::ViewKey yPartitionsOverride      = {"yPartitionsOverride"};      ///< Override of number of
                                                                                     ///< subdivisions in y key
    dataRepository::ViewKey zPartitionsOverride      = {"zPartitionsOverride"};      ///< Override of number of
                                                                                     ///< subdivisions in z key
    dataRepository::ViewKey overridePartitionNumbers = {"overridePartitionNumbers"}; ///< Flag to override partitioning
                                                                                     ///< key
    dataRepository::ViewKey schemaFileName           = {"schemaFileName"};           ///< Schema file name key
    dataRepository::ViewKey problemName              = {"problemName"};              ///< Problem name key
    dataRepository::ViewKey outputDirectory          = {"outputDirectory"};          ///< Output directory key
    dataRepository::ViewKey useNonblockingMPI        = {"useNonblockingMPI"};        ///< Flag to use non-block MPI key
    dataRepository::ViewKey suppressPinned           = {"suppressPinned"};           ///< Flag to suppress use of pinned
                                                                                     ///< memory key
  } viewKeys; ///< Command line input viewKeys

  /// Child group viewKeys
  struct groupKeysStruct
  {
    static constexpr auto numericalMethodsManagerString = "NumericalMethods";             ///< Numerical methods string
    dataRepository::GroupKey commandLine    = { "commandLine" };                          ///< Command line key
    dataRepository::GroupKey constitutiveManager = { "Constitutive" };                    ///< Constitutive key
    dataRepository::GroupKey domain    = { "domain" };                                    ///< Domain key
    dataRepository::GroupKey eventManager = { "Events" };                                 ///< Events key
    dataRepository::GroupKey fieldSpecificationManager = { "FieldSpecifications" };       ///< Field specification key
    dataRepository::GroupKey functionManager = { "Functions" };                           ///< Functions key
    dataRepository::GroupKey geometricObjectManager = { "Geometry" };                     ///< Geometry key
    dataRepository::GroupKey meshManager = { "Mesh" };                                    ///< Mesh key
    dataRepository::GroupKey numericalMethodsManager = { numericalMethodsManagerString }; ///< Numerical methods key
    dataRepository::GroupKey outputManager = { "Outputs" };                               ///< Outputs key
    dataRepository::GroupKey physicsSolverManager = { "Solvers" };                        ///< Solvers key
    dataRepository::GroupKey graphManager = { "Graphs" };                                 ///< Graphs key

  } groupKeys; ///< Child group viewKeys

  /**
   * @brief Returns the PhysicsSolverManager
   * @return Reference to the PhysicsSolverManager
   */
  PhysicsSolverManager & GetPhysicsSolverManager()
  {
    return *m_physicsSolverManager;
  }

  /**
   * @brief Returns the PhysicsSolverManager
   * @return Const reference to the PhysicsSolverManager
   */
  PhysicsSolverManager const & GetPhysicsSolverManager() const
  {
    return *m_physicsSolverManager;
  }

protected:
  /**
   * @brief Post process the command line input
   */
  virtual void PostProcessInput() override final;

private:

  /**
   * @brief Determine the number of quadrature points required for each
   *   subregion.
   * @param meshBodies Reference to the mesh bodies object.
   * @return A map containing the number of quadrature points for every
   *   region/subregion key pair.
   *
   * Checks all physics solvers for targetRegions and constitutive models to
   * determine the minimum number of quadrature points for each subregion.
   */
  map< std::pair< string, string >, localIndex > calculateRegionQuadrature( Group & meshBodies );

  /**
   * @brief Allocate constitutive relations on each subregion with appropriate
   *   number of quadrature point.
   * @param meshBodies Reference to the mesh bodies object.
   * @param constitutiveManager The constitutive manager object.
   * @param regionQuadrature The map containing the number of quadrature points for every subregion.
   */
  void setRegionQuadrature( Group & meshBodies,
                            constitutive::ConstitutiveManager const & constitutiveManager,
                            map< std::pair< string, string >, localIndex > const & regionQuadrature );

  /// The PhysicsSolverManager
  PhysicsSolverManager * m_physicsSolverManager;

  /// The EventManager
  EventManager * m_eventManager;

  /// The FunctionManager
  FunctionManager * m_functionManager;

 /// The GraphManager
  GraphManager * m_graphManager;
};

} /* namespace geosx */

#endif /* GEOSX_MANAGERS_PROBLEMMANAGER_HPP_ */
