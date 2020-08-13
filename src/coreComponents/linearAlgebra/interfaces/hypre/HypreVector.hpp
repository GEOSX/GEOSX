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
 * @file HypreVector.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_HYPREVECTOR_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_HYPREVECTOR_HPP_

#include "linearAlgebra/interfaces/VectorBase.hpp"

/**
 * @name Hypre forward declarations.
 *
 * Forward declare hypre's vector structs and pointer aliases in order
 * to avoid including hypre headers and leaking into the rest of GEOSX.
 */
///@{

/// IJVector struct forward declaration
extern "C" struct hypre_IJVector_struct;

/// IJVector pointer alias
using HYPRE_IJVector = hypre_IJVector_struct *;

/// ParVector struct forward definition
extern "C" struct hypre_ParVector_struct;

/// ParVector pointer alias
using HYPRE_ParVector = hypre_ParVector_struct *;

///@}

namespace geosx
{
/**
 * @brief Wrapper class for hypre's ParVector.
 *
 * This class creates and provides basic support for the HYPRE_ParVector object
 * type used in Hypre using the linear-algebraic system interface (IJ interface).
 */
class HypreVector final : private VectorBase< HypreVector >
{
public:
  /**
   * @name Constructor/Destructor Methods
   */
  ///@{

  /**
   * @brief Empty vector constructor.
   * Create an empty (distributed) vector.
   */
  HypreVector();

  /**
   * @brief Copy constructor.
   * @param src vector to be copied
   */
  HypreVector( HypreVector const & src );

  /**
   * @brief Move constructor.
   * @param src vector to be moved
   */
  HypreVector( HypreVector && src ) noexcept;

  /**
   * @brief Copy assignment.
   * @param src HypreVector to be copied.
   * @return the new vector.
   */
  HypreVector &
  operator=( HypreVector const & src );

  /**
   * @brief Move assignment.
   * @param src HypreVector to be moved.
   * @return the new vector.
   */
  HypreVector &
  operator=( HypreVector && src ) noexcept;

  /**
   * @brief Destructor.
   */
  ~HypreVector();

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
  set( globalIndex const globalRowIndex, real64 const value ) override;

  virtual void
  add( globalIndex const globalRowIndex, real64 const value ) override;

  virtual void
  set( globalIndex const * globalRowIndices,
       real64 const * values,
       localIndex size ) override;

  virtual void
  add( globalIndex const * globalRowIndices,
       real64 const * values,
       localIndex const size ) override;

  virtual void
  set( arraySlice1d< globalIndex const > const & globalRowIndices,
       arraySlice1d< real64 const > const & values ) override;

  virtual void
  add( arraySlice1d< globalIndex const > const & globalRowIndices,
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
  dot( HypreVector const & vec ) const override;

  virtual void
  copy( HypreVector const & x ) override;

  virtual void
  axpy( real64 const alpha, HypreVector const & x ) override;

  virtual void
  axpby( real64 const alpha,
         HypreVector const & x,
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
  get( globalIndex globalRow ) const override;

  virtual void
  get( arraySlice1d< globalIndex const > const & globalRowIndices,
       arraySlice1d< real64 > const & values ) const override;

  virtual localIndex
  getLocalRowID( globalIndex const globalRowIndex ) const override;

  virtual globalIndex
  getGlobalRowID( localIndex const localRowIndex ) const override;

  virtual real64 const *
  extractLocalVector() const override;

  virtual real64 *
  extractLocalVector() override;

  virtual MPI_Comm
  getComm() const override;

  virtual void
  print( std::ostream & os = std::cout ) const override;

  virtual void
  write( string const & filename,
         LAIOutputFormat const format =
           LAIOutputFormat::MATRIX_MARKET ) const override;

  ///@}

  /**
   * @brief Returns a pointer to the implementation.
   * @return the underlying HYPRE_ParVector object.
   */
  HYPRE_ParVector const &
  unwrapped() const;

  /**
   * @brief Returns a pointer to the implementation.
   * @return the underlying HYPRE_IJVector object.
   */
  HYPRE_IJVector const &
  unwrappedIJ() const;

private:
  /**
   * Pointer to underlying HYPRE_IJVector type.
   */
  HYPRE_IJVector m_ij_vector;

  /**
   * Pointer to underlying HYPRE_ParVector type.
   */
  HYPRE_ParVector m_par_vector;
};

}  // end namespace geosx

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_HYPREVECTOR_HPP_*/
