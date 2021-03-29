#ifndef __PROJECTIONEDFMHELPER_H_
#define __PROJECTIONEDFMHELPER_H_

#include "finiteVolume/FluxApproximationBase.hpp"  // FIXME: include appropriate explicitly
#include <list>                                  // provides std::vector

namespace geosx {

class ProjectionEDFMHelper {
 public:
  ProjectionEDFMHelper( MeshLevel const & mesh,
                        GeometricObjectManager const & geometricObjManager,
                        std::string const & coeffName,
                        CellElementStencilTPFA & stencil );

  // add Fracture-matrix connections to the cell stencil
  void addNonNeighboringConnections(EmbeddedSurfaceSubRegion const & fractureSubRegion) const;

  virtual ~ProjectionEDFMHelper() = default;

 private:
  // select cell faces that will host non-neighboring fracture-matrix connections
  std::list<localIndex> selectFaces(FixedOneToManyRelation const & subRegionFaces,
                                      CellDescriptor const & hostCellID,
                                      localIndex const fracElement,
                                      EmbeddedSurfaceSubRegion const & fractureSubRegion) const;

  // check the intersection  a fracture element and an edge
  bool intersection( real64 (&fracOrigin)[3],
                     arraySlice1d< real64 const > const & fracNormal,
                     localIndex edgeIdx,
                     real64 (&tmp)[3] ) const noexcept;

  // returns true is the face has only one neighbor
  bool isBoundaryFace( localIndex faceIdx ) const noexcept;

  // check if the center of a face is on the same side of the fracture as cell center
  bool onLargerSide( localIndex faceIdx,
                     real64 signedDistanceCellCenterToFrac,
                     real64 (&fracOrigin)[3],
                     arraySlice1d< real64 const > const & fracNormal ) const noexcept;

  // compute the signed distance between the fracture and as cell center
  real64 getSignedDistanceCellCenterToFracPlane( CellDescriptor const & hostCellID,
                                                 arraySlice1d< real64 const > const & fracNormal,
                                                 real64 const (&fracOrigin)[3],
                                                 real64 (&tmp)[3] ) const noexcept;

  // returns true if the signed distance from the neighbor center to the frac is of the same
  // sign as signedDistanceCellCenterToFrac (computed in the host cell)
  bool neighborOnSameSide( localIndex faceIdx,
                           real64 signedDistanceCellCenterToFrac,
                           CellDescriptor const & hostCellID,
                           EmbeddedSurfaceSubRegion const & fractureSubRegion ) const;

  // given a face and its neighboring cell, return the id of the other cell
  CellDescriptor otherCell( localIndex faceIdx, CellDescriptor const & hostCellID ) const;

  // compute the absolute transmissibility for non-neighboring F-M connection
  real64 fractureMatrixTransmissilibility( CellDescriptor const & neighborCell,
                                           localIndex fracElement,
                                           EmbeddedSurfaceSubRegion const & fractureSubRegion,
                                           localIndex faceIdx ) const;

  // add non-neighboring F-M connection with given transmissibility tothe cell stencil
  void addNonNeighboringConnection( localIndex fracElement,
                                    CellDescriptor const & cell,
                                    real64 transmissibility,
                                    EmbeddedSurfaceSubRegion const & fractureSubRegion ) const;

  // Private variables
  MeshLevel const & m_mesh;
  GeometricObjectManager const & m_geometricObjManager;
  ElementRegionManager const & m_elementManager;
  FaceManager const & m_faceManager;
  NodeManager const & m_nodeManager;
  EdgeManager const & m_edgeManager;
  // ArrayOfArraysView< localIndex const > const & m_faceToEdges;
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_nodesCoord;
  arrayView2d< localIndex const > const m_edgeToNodes;
  arrayView2d< localIndex const > const m_facesToCells;
  arrayView2d< localIndex const > const m_facesToRegions;
  arrayView2d< localIndex const > const m_facesToSubRegions;
  ArrayOfArraysView< localIndex const > const m_facesToNodes;
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > m_nodeReferencePosition;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > const m_cellCenters;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > const m_permTensor;
  CellElementStencilTPFA & m_stencil;
};

}  // end namespace geosx


#endif // __PROJECTIONEDFMHELPER_H_
