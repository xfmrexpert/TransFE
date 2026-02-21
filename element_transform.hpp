/***************************************************************************
 *   Copyright (C) 2005-2025 by T. C. Raymond                              *
 *   tcraymond@inductivereasoning.com                                      *
 *                                                                         *
 *   Use of this source code is governed by an MIT-style                   *
 *   license that can be found in the LICENSE.txt file or at               *
 *   https://opensource.org/licenses/MIT.                                  *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                  *
 *                                                                         *
 ***************************************************************************/

#pragma once

#include <memory>
#include "typedefs.h"
#include "Mesh/mesh.h"
#include "geometry_lagrange_shape_function.hpp"

// Forward declaration if GeometryLagrangeShapeFunction is not fully defined
class GeometryLagrangeShapeFunction;

namespace TFEM
{
    class ElementTransform {

    protected:
        std::unique_ptr<GeometryLagrangeShapeFunction> shapeFunction;
        size_t NumSpatialDims;    // Number of spatial dimensions (e.g., 2 for 2D)
        size_t NumParametricDims; // Number of parametric dimensions (e.g., 2 for quadrilateral)

    public:
        ElementTransform(size_t nsd, size_t npd, int order = 1)
            : shapeFunction(std::make_unique<GeometryLagrangeShapeFunction>(nsd, order)), NumSpatialDims(nsd), NumParametricDims(npd)
        {
        }

        virtual ~ElementTransform() = default;

        size_t referenceDimensions() const {
            return NumParametricDims;
        };

        size_t spatialDimensions() const {
            return NumSpatialDims;
        };

        virtual Matrix<double> Jacobian(const Point& pt_ref, const CellView& cell) const {
            // dG_ds: Gradients of shape function used for geometry transformation
            // [dX_dxi, dX_deta]
            // [dY_dxi, dY_deta]
            Matrix<double> dG_ds = shapeFunction->grad_N(pt_ref);

            // Initialize a nsd x npd Jacobian matrix with zeros
            Matrix<double> dXds = Matrix<double>::Zero(NumSpatialDims, NumParametricDims);

            auto nodes = cell.Nodes(); // Assumes nodes are ordered
            size_t nen = nodes.size(); // Number of element nodes

            for (size_t i = 0; i < nen; i++) { // Loop over shape functions (one per node)
                for (int j = 0; j < NumSpatialDims; j++) { // Loop over global dimension (x, y)
                    for (int k = 0; k < NumParametricDims; k++) { // Loop over local (shape) dimension (xi, eta)
                        dXds(j, k) += dG_ds(i, k) * nodes[i][j];
                    }
                }
            }

            return dXds;

        };

        virtual Point mapReferencePointToPhysical(const Point& ptRef, const CellView& cell) const
        {
            Point pt;
            const auto& N_ref = shapeFunction->N(ptRef);
            std::vector<Point> cell_nodes;
            cell.getNodes(cell_nodes);
            for (size_t i = 0; i < cell_nodes.size(); i++) {
                pt += N_ref(i) * cell_nodes[i];
            }
            return pt;
        };
    };
}

