

#ifndef GEOSX_FINITEELEMENT_FINITEELEMENTBASE_HPP_
#define GEOSX_FINITEELEMENT_FINITEELEMENTBASE_HPP_

#include "common/DataTypes.hpp"
#include "common/GeosxMacros.hpp"

namespace geosx
{
namespace finiteElement
{
/**
 * @class FiniteElementShapeFunctionKernelBase
 * @brief Base class for the finite element kernels.
 */
class FiniteElementBase
{
public:
  /**
   * @brief Virtual getter for the number of quadrature points per element.
   * @return The number of quadrature points per element.
   */
  virtual localIndex
  getNumQuadraturePoints() const = 0;

  /**
   * @brief Virtual getter for the number of support points per element.
   * @return The number of support points per element.
   */
  virtual localIndex
  getNumSupportPoints() const = 0;

  /**
   * @brief Destructor
   */
  virtual ~FiniteElementBase() = default;

  /**
   * @brief Computes the inverse of a 3x3 c-array.
   * @param J The array to invert...which is also used to store the inverse.
   * @return The determinant of @p J.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static real64 inverse( real64 ( &J )[3][3] )
  {
    real64 const temp[3][3] = { { J[1][1] * J[2][2] - J[1][2] * J[2][1],
                                  J[0][2] * J[2][1] - J[0][1] * J[2][2],
                                  J[0][1] * J[1][2] - J[0][2] * J[1][1] },
                                { J[1][2] * J[2][0] - J[1][0] * J[2][2],
                                  J[0][0] * J[2][2] - J[0][2] * J[2][0],
                                  J[0][2] * J[1][0] - J[0][0] * J[1][2] },
                                { J[1][0] * J[2][1] - J[1][1] * J[2][0],
                                  J[0][1] * J[2][0] - J[0][0] * J[2][1],
                                  J[0][0] * J[1][1] - J[0][1] * J[1][0] } };

    real64 const det =
      J[0][0] * temp[0][0] + J[1][0] * temp[0][1] + J[2][0] * temp[0][2];
    real64 const invDet = 1.0 / det;

    for( int i = 0; i < 3; ++i )
    {
      for( int j = 0; j < 3; ++j )
      {
        J[i][j] = temp[i][j] * invDet;
      }
    }
    return det;
  }

  /**
   * @brief Calculate the determinant of a 3x3 c-array.
   * @param J The input array.
   * @return The determinant of @p J
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static real64
  detJ( real64 const ( &J )[3][3] )
  {
    return J[0][0] * ( J[1][1] * J[2][2] - J[1][2] * J[2][1] ) +
      J[1][0] * ( J[0][2] * J[2][1] - J[0][1] * J[2][2] ) +
      J[2][0] * ( J[0][1] * J[1][2] - J[0][2] * J[1][1] );
  }

  //TODO we want to keep views and provide interfaces to this data here for cases
  //     where we pre-compute the shape function derivatives...maybe...tbd.
  //private:
  //  arrayView4d< real64 const > const m_dNdX;
  //  arrayView2d< real64 const > const m_detJ;
};

}  // namespace finiteElement
}  // namespace geosx

#endif  //GEOSX_CORE_FINITEELEMENT_FINITEELEMENTBASE
