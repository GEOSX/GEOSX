/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ProblemManager.hpp
 */


#ifndef GEOSX_MANAGERS_PROBLEMMANAGER_HPP_
#define GEOSX_MANAGERS_PROBLEMMANAGER_HPP_

#include "EventManager.hpp"
#include "fileIO/schema/schemaUtilities.hpp"

namespace geosx
{

class PhysicsSolverManager;
class DomainPartition;
class GeometricObjectManager;
namespace constitutive
{
class ConstitutiveManager;
}
class FunctionManager;
class FieldSpecificationManager;
struct CommandLineOptions;

/**
 * @class ProblemManager
 * @brief This is the class handling the operation flow of the problem being ran in GEOSX
 */
class ProblemManager : public dataRepository::Group
{
public:

  /**
   * @brief Create a new ProblemManager, it must be created from the root conduit node.
   * @param root The root conduit node.
   */
  explicit ProblemManager( conduit::Node & root );

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
  virtual void setSchemaDeviations( xmlWrapper::xmlNode schemaRoot,
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
  virtual Group * createChild( string const & childKey, string const & childName ) override;

  /**
   * @brief Parses command line input
   */
  void parseCommandLineInput();

  /**
   * @brief Parses a restart file.
   * @param restartFileName The name of the restart file.
   * @param options The command line options.
   * @return Flag indicating beginFromRestart status
   */
  static bool parseRestart( string & restartFileName, CommandLineOptions const & options );

  /**
   * @brief Generates the xml schema documentation
   * This function is called when the code is called with the -s schema_name option.
   * @details Before generating the schema, the code builds up a comprehensive datastructure.
   * (Note: catalog objects throughout the code will typically be registered via the
   * ExpandObjectCatalogs method.)  Once ready, SchemaUtilities will recusively walk
   * through the database, generating the xml schema.
   */
  void generateDocumentation();

  /**
   * @brief Parses the input xml file
   * @details The name of the input file is indicated via the -i option on the command line
   */
  void parseInputFile();

  /**
   * @brief Generates numerical meshes used throughout the code
   */
  void generateMesh();

  /**
   * @brief Allocates constitutive relations according to the discretizations
   *   on each subregion.
   */
  void applyNumericalMethods();

  /**
   * @brief Defines the order in which objects should be initialized
   * @param order list defining ordering sequence
   */
  void initializationOrder( string_array & order ) override final;

  /**
   * @brief Sets up the problem after the input has been read in
   */
  void problemSetup();

  /**
   * @brief Run the events in the scheduler.
   * @return True iff the simulation exited early, and needs to be run again to completion.
   */
  bool runSimulation();

  /**
   * @brief After initialization, overwrites data using a restart file
   */
  void readRestartOverwrite();

  /**
   * @brief Applies initial conditions indicated within the input file FieldSpecifications block
   */
  void applyInitialConditions();

  /**
   * @brief Returns a pointer to the DomainPartition
   * @return Pointer to the DomainPartition
   */
  DomainPartition & getDomainPartition();

  /**
   * @brief Returns a pointer to the DomainPartition
   * @return Const pointer to the DomainPartition
   */
  DomainPartition const & getDomainPartition() const;

  /**
   * @brief Returns a pointer to the GeometricObjectManager
   * @return Const reference to the GeometricObjectManager
   */
  GeometricObjectManager const & getGeometricObjectManager() const;

  /**
   * @brief Returns the problem name
   * @return The problem name
   */
  string const & getProblemName() const
  { return getGroup< Group >( groupKeys.commandLine ).getReference< string >( viewKeys.problemName ); }

  /**
   * @brief Returns the input file name
   * @return The input file name
   */
  string const & getInputFileName() const
  { return getGroup< Group >( groupKeys.commandLine ).getReference< string >( viewKeys.inputFileName ); }

  /**
   * @brief Returns the restart file name
   * @return The restart file name
   */
  string const & getRestartFileName() const
  { return getGroup< Group >( groupKeys.commandLine ).getReference< string >( viewKeys.restartFileName ); }

  /**
   * @brief Returns the schema file name
   * @return The schema file name
   */
  string const & getSchemaFileName() const
  { return getGroup< Group >( groupKeys.commandLine ).getReference< string >( viewKeys.schemaFileName ); }

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
    /// @return Numerical methods string
    static constexpr char const * numericalMethodsManagerString() { return "NumericalMethods"; }
    dataRepository::GroupKey commandLine    = { "commandLine" };                          ///< Command line key
    dataRepository::GroupKey constitutiveManager = { "Constitutive" };                    ///< Constitutive key
    dataRepository::GroupKey domain    = { "domain" };                                    ///< Domain key
    dataRepository::GroupKey eventManager = { "Events" };                                 ///< Events key
    dataRepository::GroupKey fieldSpecificationManager = { "FieldSpecifications" };       ///< Field specification key
    dataRepository::GroupKey functionManager = { "Functions" };                           ///< Functions key
    dataRepository::GroupKey geometricObjectManager = { "Geometry" };                     ///< Geometry key
    dataRepository::GroupKey meshManager = { "Mesh" };                                    ///< Mesh key
    dataRepository::GroupKey numericalMethodsManager = { numericalMethodsManagerString() }; ///< Numerical methods key
    dataRepository::GroupKey outputManager = { "Outputs" };                               ///< Outputs key
    dataRepository::GroupKey physicsSolverManager = { "Solvers" };                        ///< Solvers key
    dataRepository::GroupKey tasksManager = { "Tasks" };                                  ///< Tasks key
  } groupKeys; ///< Child group viewKeys

  /**
   * @brief Returns the PhysicsSolverManager
   * @return Reference to the PhysicsSolverManager
   */
  PhysicsSolverManager & getPhysicsSolverManager()
  {
    return *m_physicsSolverManager;
  }

  /**
   * @brief Returns the PhysicsSolverManager
   * @return Const reference to the PhysicsSolverManager
   */
  PhysicsSolverManager const & getPhysicsSolverManager() const
  {
    return *m_physicsSolverManager;
  }

  /**
   * @brief Returns the FunctionManager.
   * @return The FunctionManager.
   */
  FunctionManager & getFunctionManager()
  {
    GEOSX_ERROR_IF( m_functionManager == nullptr, "Not initialized." );
    return *m_functionManager;
  }

  /**
   * @brief Returns the const FunctionManager.
   * @return The const FunctionManager.
   */
  FunctionManager const & getFunctionManager() const
  {
    GEOSX_ERROR_IF( m_functionManager == nullptr, "Not initialized." );
    return *m_functionManager;
  }

  /**
   * @brief Returns the FieldSpecificationManager.
   * @return The FieldSpecificationManager.
   */
  FieldSpecificationManager & getFieldSpecificationManager()
  {
    GEOSX_ERROR_IF( m_fieldSpecificationManager == nullptr, "Not initialized." );
    return *m_fieldSpecificationManager;
  }

  /**
   * @brief Returns the const FunctionManager.
   * @return The const FunctionManager.
   */
  FieldSpecificationManager const & getFieldSpecificationManager() const
  {
    GEOSX_ERROR_IF( m_fieldSpecificationManager == nullptr, "Not initialized." );
    return *m_fieldSpecificationManager;
  }

protected:
  /**
   * @brief Post process the command line input
   */
  virtual void postProcessInput() override final;

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

  /// The FieldSpecificationManager
  FieldSpecificationManager * m_fieldSpecificationManager;
};

} /* namespace geosx */

#endif /* GEOSX_MANAGERS_PROBLEMMANAGER_HPP_ */
