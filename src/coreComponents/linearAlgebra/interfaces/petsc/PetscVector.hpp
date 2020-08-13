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
 * @file PetscVector.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_PETSCVECTOR_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_PETSCVECTOR_HPP_

#include "linearAlgebra/interfaces/VectorBase.hpp"

/**
 * @name PETSc forward declarations.
 *
 * Forward declare PETSc's vector struct and pointer aliases in order
 * to avoid including PETSc headers and leaking into the rest of GEOSX.
 */
///@{

/// Vec struct forward declaration
extern "C" struct _p_Vec;

///@}

namespace geosx
{
/**
 * @brief This class creates and provides basic support for Vec
 *        vector object type used in PETSc.
 */
class PetscVector final : private VectorBase< PetscVector >
{
public:
  /// Alias for PETSc vector struct pointer
  using Vec = struct _p_Vec *;

  /**
   * @name Constructor/Destructor Methods
   */
  ///@{

  /**
   * @brief Empty vector constructor.
   */
  PetscVector();

  /**
   * @brief Copy constructor.
   * @param src PetscVector to be copied.
   */
  PetscVector( PetscVector const & src );

  /**
   * @brief Move constructor
   * @param src PetscVector to move from
   */
  PetscVector( PetscVector && src ) noexcept;

  /**
   * @brief Copy assignment.
   * @param src PetscVector to be copied.
   * @return the new vector.
   */
  PetscVector &
  operator=( PetscVector const & src );

  /**
   * @brief Move assignment.
   * @param src PetscVector to be moved from.
   * @return the new vector.
   */
  PetscVector &
  operator=( PetscVector && src ) noexcept;

  /**
   * @brief Destructor.
   */
  ~PetscVector();

  ///@}

  /**
   * @name VectorBase interface
   */
  ///@{

  using VectorBase::closed;
  using VectorBase::extract;
  using VectorBase::ready;

  virtual bool
  created() const override;

  virtual void
  createWithLocalSize( localIndex const localSize,
                       MPI_Comm const & comm ) override;

  virtual void
  createWithGlobalSize( globalIndex const globalSize,
                        MPI_Comm const & comm ) override;

  virtual void
  create( arrayView1d< real64 const > const & localValues,
          MPI_Comm const & comm ) override;

  virtual void
  open() override;

  virtual void
  close() override;

  virtual void
  reset() override;

  virtual void
  set( globalIndex const globalRow, real64 const value ) override;

  virtual void
  add( globalIndex const globalRow, real64 const value ) override;

  virtual void
  set( globalIndex const * globalIndices,
       real64 const * values,
       localIndex size ) override;

  virtual void
  add( globalIndex const * globalIndices,
       real64 const * values,
       localIndex size ) override;

  virtual void
  set( arraySlice1d< globalIndex const > const & globalIndices,
       arraySlice1d< real64 const > const & values ) override;

  virtual void
  add( arraySlice1d< globalIndex const > const & globalIndices,
       arraySlice1d< real64 const > const & values ) override;

  virtual void
  set( real64 const value ) override;

  virtual void
  zero() override;

  virtual void
  rand( unsigned const seed = 1984 ) override;

  virtual void
  scale( real64 const scalingFactor ) override;

  virtual void
  reciprocal() override;

  virtual real64
  dot( PetscVector const & vec ) const override;

  virtual void
  copy( PetscVector const & x ) override;

  virtual void
  axpy( real64 const alpha, PetscVector const & x ) override;

  virtual void
  axpby( real64 const alpha,
         PetscVector const & x,
         real64 const beta ) override;

  virtual real64
  norm1() const override;

  virtual real64
  norm2() const override;

  virtual real64
  normInf() const override;

  virtual globalIndex
  globalSize() const override;

  virtual localIndex
  localSize() const override;

  virtual globalIndex
  ilower() const override;

  virtual globalIndex
  iupper() const override;

  virtual real64
  get( globalIndex const globalRow ) const override;

  void
  get( arraySlice1d< globalIndex const > const & globalIndices,
       arraySlice1d< real64 > const & values ) const override;

  virtual MPI_Comm
  getComm() const override;

  virtual void
  print( std::ostream & os = std::cout ) const override;

  virtual void
  write( string const & filename,
         LAIOutputFormat const format =
           LAIOutputFormat::MATRIX_MARKET ) const override;

  virtual localIndex
  getLocalRowID( globalIndex const globalRow ) const override;

  virtual globalIndex
  getGlobalRowID( localIndex const localRow ) const override;

  virtual real64 const *
  extractLocalVector() const override;

  virtual real64 *
  extractLocalVector() override;

  ///@}

  /**
   * @brief Returns a const pointer to the underlying Vec.
   * @return the const pointer to the underlying Vec.
   */
  const Vec &
  unwrapped() const;

  /**
   * @brief Returns a non-const pointer to the underlying Vec.
   * @return the non-const pointer to the underlying Vec.
   */
  Vec &
  unwrapped();

protected:
  /**
   * Pointer to underlying PETSc Vec
   */
  Vec m_vec;
};

}  // namespace geosx

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_PETSCVECTOR_HPP_*/
