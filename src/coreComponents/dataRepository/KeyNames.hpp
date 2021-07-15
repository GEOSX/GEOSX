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
 * @file KeyNames.hpp
 */

#ifndef GEOSX_DATAREPOSITORY__KEYNAMES_HPP_
#define GEOSX_DATAREPOSITORY__KEYNAMES_HPP_

#include <string>

namespace geosx
{
namespace dataRepository
{
namespace keys
{

/// @cond DO_NOT_DOCUMENT

//static constexpr auto ReferencePosition = "ReferencePosition";
static constexpr auto referencePositionString = "ReferencePosition";
static constexpr auto Temperature = "Temperature";
static constexpr auto IncrementalTemperature = "IncrementalTemperature";
static constexpr auto TotalDisplacement = "TotalDisplacement";
static constexpr auto IncrementalDisplacement = "IncrementalDisplacement";
static constexpr auto Velocity = "Velocity";
static constexpr auto Acceleration = "Acceleration";
static constexpr auto Mass = "Mass";
static constexpr auto Force = "Force";
static constexpr auto Strain = "Strain";
static constexpr auto Name = "name";
static constexpr auto Size = "size";
static constexpr auto ProblemManager = "Problem";
static constexpr auto ConstitutiveManager = "Constitutive";
static constexpr auto ConstitutiveBase = "ConstitutiveBase";
static constexpr auto solverNames = "solverNames";

static constexpr auto schema = "schema";

static constexpr auto time = "time";
static constexpr auto cycle = "cycle";
static constexpr auto beginTime = "beginTime";
static constexpr auto endTime = "endTime";
static constexpr auto dt = "dt";

static constexpr auto domain  = "domain";
static constexpr auto solvers = "solvers";
static constexpr auto simulationParameterMap = "simulationParameterMap";
static constexpr auto FE_Space    = "FE_Space";
//static constexpr auto FEM_Nodes    = "FEM_Nodes";
//static constexpr auto FEM_Edges    = "FEM_Edges";
//static constexpr auto FEM_Faces    = "FEM_Faces";
//static constexpr auto FEM_Elements = "FEM_Elements";
static constexpr auto cellManager = "cellManager";
static constexpr auto functionManager = "FunctionManager";

/// @endcond

}
}
}
#endif /* GEOSX_DATAREPOSITORY__KEYNAMES_HPP_ */
