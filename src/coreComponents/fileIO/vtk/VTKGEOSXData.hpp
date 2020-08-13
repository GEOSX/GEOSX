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

#ifndef GEOSX_FILEIO_VTK_VTKGEOSXDATA_HPP_
#define GEOSX_FILEIO_VTK_VTKGEOSXDATA_HPP_

#include "common/DataTypes.hpp"

#include <vtkAOSDataArrayTemplate.h>

namespace geosx
{
namespace vtk
{
/*!
 * @brief VTK GEOSX data class.
 * @details This class let us to deal with special types such
 * as R1Tensor to be output.
 */
class VTKGEOSXData : public vtkAOSDataArrayTemplate< real64 >
{
public:
  /*!
   * @brief Factory function
   * @return VTK GEOSX data class
   */
  static VTKGEOSXData *
  New()
  {
    VTK_STANDARD_NEW_BODY( VTKGEOSXData );
  }
  /*!
   * @brief insert Value in a custom way for R1Tensors
   * @param[in] index index where the value \p val will be inserted
   * @param[in] val value to be inserted
   */
  template< typename T >
  void
  CustomInsertValue( localIndex index, T const & val )
  {
    this->InsertValue( index, val );
  }
};

/**
 * @brief Custom insert function for an R1Tensor.
 * @param[in] index position index where the value \p val will be inserted
 * @param[in] val R1Tensor to be inserted
 */
template<>
void
VTKGEOSXData::CustomInsertValue< R1Tensor >( localIndex index,
                                             R1Tensor const & val );

}  // namespace vtk
}  // namespace geosx

#endif
