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
 * @file UniformStrainHexahedron.h
 */

#include "legacy/ElementLibrary/FiniteElement.h"

#ifndef UNIFORMSTRAINHEXAHEDRON_H_
  #define UNIFORMSTRAINHEXAHEDRON_H_

class UniformStrainHexahedron : public FiniteElement< 3 >
{
public:
  UniformStrainHexahedron();
  virtual ~UniformStrainHexahedron();

  void
  reinit( const std::vector< R1TensorT< 3 > > & mapped_support_points );

  void
  zero_energy_mode_control( const array1d< R1Tensor > & dNdx,
                            const realT & volume,
                            const array1d< R1Tensor > & x,
                            const array1d< R1Tensor > & vel,
                            const realT & dampcoef,
                            const realT & stiffcoef,
                            const realT & rho,
                            const realT & modulus,
                            const realT & dt,
                            array1d< R1Tensor > & Qstiffness,
                            array1d< R1Tensor > & force );
};

#endif /* UNIFORMSTRAINHEXAHEDRON_H_ */
