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
 * @file EpetraVector.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_EPETRAVECTOR_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_EPETRAVECTOR_HPP_

#include "linearAlgebra/interfaces/VectorBase.hpp"

class Epetra_FEVector;

namespace geosx
{
/**
 * @brief This class creates and provides basic support for the Epetra_FEVector
 *        vector object type used in Trilinos.  We use the FE version because
 *        Epetra_Vector support for long globalIDs is haphazard.
 */
class EpetraVector final : private VectorBase< EpetraVector >
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
  EpetraVector();

  /**
   * @brief Copy constructor.
   * @param src EpetraVector to be copied.
   */
  EpetraVector( EpetraVector const & src );

  /**
   * @brief Move constructor
   * @param src EpetraVector to move from
   */
  EpetraVector( EpetraVector && src ) noexcept;

  /**
   * @brief Copy assignment.
   * @param src EpetraVector to be copied.
   * @return the new vector
   */
  EpetraVector &
  operator=( EpetraVector const & src );

  /**
   * @brief Move assignment.
   * @param src EpetraVector to be moved from.
   * @return the new vector
   */
  EpetraVector &
  operator=( EpetraVector && src ) noexcept;

  /**
   * @brief Destructor.
   */
  ~EpetraVector();

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
  dot( EpetraVector const & vec ) const override;

  virtual void
  copy( EpetraVector const & x ) override;

  virtual void
  axpy( real64 const alpha, EpetraVector const & x ) override;

  virtual void
  axpby( real64 const alpha,
         EpetraVector const & x,
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
  get( arraySlice1d< globalIndex const > const & globalIndices,
       arraySlice1d< real64 > const & values ) const override;

  virtual localIndex
  getLocalRowID( globalIndex const globalRow ) const override;

  virtual globalIndex
  getGlobalRowID( localIndex const localRow ) const override;

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
   * @brief Returns a const pointer to the underlying Epetra object.
   * @return const pointer to the underlying Epetra object
   */
  Epetra_FEVector const &
  unwrapped() const;

  /**
   * @brief Returns a non-const pointer to the underlying Epetra object.
   * @return non-const pointer to the underlying Epetra object
   */
  Epetra_FEVector &
  unwrapped();

private:
  /// Unique pointer to underlying Epetra_FEVector object.
  std::unique_ptr< Epetra_FEVector > m_vector;
};

}  // namespace geosx

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_EPETRAVECTOR_HPP_*/
