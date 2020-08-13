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
 * @file HypreMGRStrategies.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_HYPREMGRSTRATEGIES_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_HYPREMGRSTRATEGIES_HPP_

#include "common/DataTypes.hpp"

#include <_hypre_utilities.h>

namespace geosx
{
/**
 * @todo Add a detailed description with an example
 *
 * @brief Compute an array of unique component labels.
 *
 * @details
 *
 * @param numComponentsPerField array of number of components per field
 * @param numLocalDofsPerField array of local number of dofs per field
 * @return array1d of HYPRE_Int labels
 */
inline array1d< HYPRE_Int >
computeLocalDofComponentLabels(
  arraySlice1d< localIndex const > const & numComponentsPerField,
  arraySlice1d< localIndex const > const & numLocalDofsPerField )
{
  array1d< HYPRE_Int > ret;
  HYPRE_Int numFields =
    LvArray::integerConversion< HYPRE_Int >( numLocalDofsPerField.size() );

  if( numFields > 0 )
  {
    HYPRE_Int numTotalLocalDof = 0;
    for( HYPRE_Int i = 0; i < numFields; ++i )
    {
      numTotalLocalDof +=
        LvArray::integerConversion< HYPRE_Int >( numLocalDofsPerField[i] );
    }

    ret.resize( numTotalLocalDof );

    HYPRE_Int firstLabel = 0;
    HYPRE_Int istr = 0;
    HYPRE_Int iend;
    HYPRE_Int numComp;
    for( HYPRE_Int iFld = 0; iFld < numFields; ++iFld )
    {
      numComp =
        LvArray::integerConversion< HYPRE_Int >( numComponentsPerField[iFld] );
      array1d< HYPRE_Int > vectorLabels( numComp );
      for( HYPRE_Int k = 0; k < numComp; ++k )
      {
        vectorLabels[k] = k + firstLabel;
      }
      iend =
        istr + LvArray::integerConversion< HYPRE_Int >( numLocalDofsPerField[iFld] );
      ;
      for( localIndex i = istr; i < iend; i += numComp )
      {
        for( integer k = 0; k < numComp; ++k )
        {
          ret[i + k] = vectorLabels[k];
        }
      }
      istr += iend;
      firstLabel += numComp;
    }
  }
  return ret;
}

}  // namespace geosx

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_HYPREMGRSTRATEGIES_HPP_*/
